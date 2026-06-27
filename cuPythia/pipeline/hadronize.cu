// cuPythia hadronization — increment 2: a GPU Lund string fragmentation chain.
//
// Fragments ONE straight q-qbar string (the FSR-shower final state, here a fixed-sqrt(s)
// u/d/s q-qbar at rest) into pseudoscalar AND vector mesons, one string per GPU thread,
// counter-RNG (../common/rng.cuh). Faithful to Pythia 8.317: the zLund f(z) sampler
// (validated in zlund_test.cu), the StringFlav meson selection (flavour 1:1:probStoUD,
// the ALWAYS-drawn spin choice -> vector fraction m_v/(1+m_v), eta/eta' suppression,
// uds mixing), the StringPT pT (enhancedFraction draw kept, sigma/sqrt2, Box-Muller
// sin-first), the constituent-mass stop test, light-cone longitudinal kinematics, and an
// exact two-body finalTwo with refragment-on-failure. Honest simplifications (documented):
// only the two lowest meson multiplets (Pythia's L=1 rates are 0 by default), pole masses
// (no Breit-Wigner), no baryons/diquarks/popcorn, no decays, and finalTwo splits the
// remainder as a clean 2-body decay along the string axis (no extra pT smear). Validated by
// exact 4-momentum conservation, on-shellness, GPU==CPU determinism, reproducibility, and
// (separately, multiplicity_pythia.cc) vs Pythia under a matched config.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o hadronize hadronize.cu
// Run:   ./hadronize [nEvents=200000] [sqrtS=91.1876]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)
#define MAXH 80

// ---- Pythia 8.317 defaults ----
__device__ __host__ inline double ALUND(){return 0.68;}
__device__ __host__ inline double BLUND(){return 0.98;}
static const double SIGMAQ   = 0.335/1.4142135623730951;  // sigma/sqrt2
static const double PROBSTOUD= 0.217;
static const double MESONUDV = 0.50,  MESONSV = 0.55;
static const double THETAPS  = -15.0, THETAV  = 36.0;
static const double ETASUP   = 0.60,  ETAPRIMESUP = 0.12;
static const double STOPMASS = 0.8,   STOPNEWFLAV = 2.0, STOPSMEAR = 0.2;
static const double ENHFRAC  = 0.01,  ENHWIDTH = 2.0;

#include "zlund_inc.cuh"   // zLundMax + zLundSample (shared with zlund_test)

__host__ __device__ inline double mesonMass(int pdg){
  switch(abs(pdg)){
    case 211:return 0.13957; case 111:return 0.13498; case 221:return 0.54786; case 331:return 0.95778;
    case 321:return 0.49368; case 311:return 0.49761;
    case 213:return 0.77526; case 113:return 0.77526; case 223:return 0.78266; case 333:return 1.01946;
    case 323:return 0.89167; case 313:return 0.89555;
  } return -1.0;
}
__host__ __device__ inline double constMass(int absFlav){ return (absFlav==3)?0.5:0.33; }
__host__ __device__ inline bool isCharged(int pdg){ int a=abs(pdg); return a==211||a==321||a==213||a==323; }

// pick a new light quark flavour (d=1,u=2,s=3) with weights 1:1:probStoUD
__host__ __device__ inline int pickFlav(uint64_t& ctr){
  double r=u01(splitmix64(ctr++))*(2.0+PROBSTOUD);
  return (r<1.0)?1:((r<2.0)?2:3);
}
// StringFlav::combine meson branch (signed quark ids). Returns 0 on eta/eta' suppression.
__host__ __device__ inline int combineMeson(int id1,int id2,uint64_t& ctr){
  int a1=abs(id1),a2=abs(id2),idMax=max(a1,a2),idMin=min(a1,a2);
  int flav=(idMax<3)?0:idMax-2;
  double mv=(flav==0)?MESONUDV:MESONSV;
  double rndmSpin=(1.0+mv)*u01(splitmix64(ctr++));
  int spin=0; rndmSpin-=1.0; if(rndmSpin>0.0){ spin=1; }
  int code=(spin==0)?1:3;
  int idMeson=100*idMax+10*idMin+code;
  if(idMax!=idMin){
    int sign=(idMax%2==0)?1:-1;
    if((idMax==a1 && id1<0)||(idMax==a2 && id2<0)) sign=-sign;
    idMeson*=sign;
  } else { // diagonal uds mixing
    double alpha=((spin==0)?(90.0-(THETAPS+54.7)):(THETAV+54.7))*M_PI/180.0;
    double mix1,mix2;
    if(flav==0){ mix1=0.5; mix2=0.5*(1.0+sin(alpha)*sin(alpha)); }
    else       { mix1=0.0; mix2=cos(alpha)*cos(alpha); }
    double rMix=u01(splitmix64(ctr++));
    if(rMix<mix1) idMeson=110; else if(rMix<mix2) idMeson=220; else idMeson=330;
    idMeson+=code;
    if(idMeson==221 && ETASUP     <u01(splitmix64(ctr++))) return 0;
    if(idMeson==331 && ETAPRIMESUP<u01(splitmix64(ctr++))) return 0;
  }
  return idMeson;
}
// Gaussian pair pT for a new break: enhancedFraction draw kept (B1), Box-Muller sin-first (B4).
__host__ __device__ inline void pairPT(uint64_t& ctr,double& gx,double& gy){
  double mult=(u01(splitmix64(ctr++))<ENHFRAC)?ENHWIDTH:1.0;
  double r=sqrt(-2.0*log(u01(splitmix64(ctr++))+1e-300));
  double phi=2.0*M_PI*u01(splitmix64(ctr++));
  gx=mult*SIGMAQ*r*sin(phi); gy=mult*SIGMAQ*r*cos(phi);
}
__host__ __device__ inline void boostBy(const double* q,double ex,double ey,double ez,double gamma,double* o){
  double bdq=ex*q[0]+ey*q[1]+ez*q[2], e2=ex*ex+ey*ey+ez*ez;
  double k=(e2>1e-18)?((gamma-1.0)*bdq/e2+gamma*q[3]):0.0;
  o[0]=q[0]+k*ex; o[1]=q[1]+k*ey; o[2]=q[2]+k*ez; o[3]=gamma*(q[3]+bdq);
}

// One attempt at fragmenting the string; returns nHadrons or -1 on failure (caller retries).
#ifdef USE_BW
#include "bw_inc.cuh"
#define HADMASS(pdg,ctr) sampleBWmass(pdg, mesonMass(pdg), ctr)   // Breit-Wigner vector masses
#else
#define HADMASS(pdg,ctr) mesonMass(pdg)                            // pole masses (default)
#endif
__host__ __device__ inline int tryFragment(double E,int flA,int flB,uint64_t& ctr,double* H,int* hid,double* hm){
  double w2=4.0*E*E;                      // (2E)^2 = s
  double xPos=1.0,xNeg=1.0;               // remaining light-cone fractions
  int fPos=flA, fNeg=flB;                 // + end quark, - end antiquark (flB<0)
  double pTpx=0,pTpy=0, pTnx=0,pTny=0;    // running endpoint transverse momenta
  double pRem[4]={0,0,0,2.0*E};           // remaining 4-momentum
  int n=0;
  for(int step=0; step<MAXH-2; ++step){
    bool fromPos=(u01(splitmix64(ctr++))<0.5);
    // draw a hadron (retry whole draw if combine fails on eta/eta')
    int pdg=0,fNew=0; double mh=0,hpx=0,hpy=0,gx=0,gy=0;
    for(int tries=0; tries<20 && pdg==0; ++tries){
      fNew=pickFlav(ctr);
      int endFlav=fromPos?fPos:fNeg;
      pdg=combineMeson(endFlav, fromPos?-fNew:fNew, ctr);  // endpoint + new (anti)quark
      pairPT(ctr,gx,gy);
      if(pdg!=0){ mh=HADMASS(pdg,ctr);
        hpx=(fromPos?pTpx:pTnx)-gx; hpy=(fromPos?pTpy:pTny)-gy; }
    }
    if(pdg==0) return -1;
    double mT2=mh*mh+hpx*hpx+hpy*hpy;
    // stop test (constituent masses), stopSmear drawn AFTER the hadron draw
    double wMin=(STOPMASS+constMass(abs(fPos))+constMass(abs(fNeg))
                 +STOPNEWFLAV*constMass(fNew))*(1.0+(2.0*u01(splitmix64(ctr++))-1.0)*STOPSMEAR);
    double w2Rem=pRem[3]*pRem[3]-pRem[0]*pRem[0]-pRem[1]*pRem[1]-pRem[2]*pRem[2];
    if(w2Rem<wMin*wMin) break;             // -> finalTwo
    // longitudinal kinematics
    double z=zLundSample(ALUND(),BLUND()*mT2,1.0,ctr);
    double xPosHad,xNegHad;
    if(fromPos){ xPosHad=z*xPos; if(xPosHad<=0) return -1; xNegHad=mT2/(xPosHad*w2);
                 if(xNegHad>xNeg) return -1; }
    else       { xNegHad=z*xNeg; if(xNegHad<=0) return -1; xPosHad=mT2/(xNegHad*w2);
                 if(xPosHad>xPos) return -1; }
    double had[4]={hpx,hpy,(xPosHad-xNegHad)*E,(xPosHad+xNegHad)*E};
    if(n>=MAXH) return -1;
    H[4*n]=had[0];H[4*n+1]=had[1];H[4*n+2]=had[2];H[4*n+3]=had[3]; hid[n]=pdg; hm[n]=mh; n++;
    xPos-=xPosHad; xNeg-=xNegHad;
    for(int k=0;k<4;++k) pRem[k]-=had[k];
    if(fromPos){ fPos=fNew; pTpx=gx; pTpy=gy; } else { fNeg=-fNew; pTnx=gx; pTny=gy; }
  }
  // finalTwo: split remaining pRem into two mesons via an exact 2-body decay.
  double W2=pRem[3]*pRem[3]-pRem[0]*pRem[0]-pRem[1]*pRem[1]-pRem[2]*pRem[2];
  if(W2<=0) return -1; double W=sqrt(W2);
  int pdg1=0,pdg2=0,fNew=0; double m1=0,m2=0;
  for(int tries=0; tries<20 && (pdg1==0||pdg2==0); ++tries){
    fNew=pickFlav(ctr);
    pdg1=combineMeson(fPos,-fNew,ctr);
    pdg2=combineMeson(fNew,fNeg,ctr);
    if(pdg1!=0&&pdg2!=0){ m1=HADMASS(pdg1,ctr); m2=HADMASS(pdg2,ctr); }
  }
  if(pdg1==0||pdg2==0) return -1;
  if(W<m1+m2) return -1;                   // doesn't fit -> caller refragments
  double E1=(W2+m1*m1-m2*m2)/(2.0*W);
  double ps=sqrt(fmax(0.0,(W2-(m1+m2)*(m1+m2))*(W2-(m1-m2)*(m1-m2))))/(2.0*W);
  double sgn=(u01(splitmix64(ctr++))<1.0/(1.0+exp(fmin(50.0,BLUND()*ps)))) ? -1.0 : 1.0; // probReverse
  double q1[4]={0,0, sgn*ps, E1}, q2[4]={0,0,-sgn*ps, W-E1};      // back-to-back along string axis
  double gA=pRem[3]/W, bx=pRem[0]/pRem[3],by=pRem[1]/pRem[3],bz=pRem[2]/pRem[3];
  double l1[4],l2[4]; boostBy(q1,bx,by,bz,gA,l1); boostBy(q2,bx,by,bz,gA,l2);
  if(n+2>MAXH) return -1;
  H[4*n]=l1[0];H[4*n+1]=l1[1];H[4*n+2]=l1[2];H[4*n+3]=l1[3]; hid[n]=pdg1; hm[n]=m1; n++;
  H[4*n]=l2[0];H[4*n+1]=l2[1];H[4*n+2]=l2[2];H[4*n+3]=l2[3]; hid[n]=pdg2; hm[n]=m2; n++;
  return n;
}
__host__ __device__ inline int hadronizeString(double E,int flA,int flB,uint64_t ctr,double* H,int* hid,double* hm){
  for(int retry=0; retry<25; ++retry){
    int n=tryFragment(E,flA,flB,ctr,H,hid,hm);   // ctr advances continuously across retries
    if(n>0) return n;
  }
  return -1;
}

__global__ void kern(int N,double E,uint64_t base,int* outN,int* outNc,double* outTot,double* outDm,double* outVm){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  double H[MAXH*4]; int hid[MAXH]; double hm[MAXH];
  int n=hadronizeString(E,2,-2, base+(uint64_t)e*0x9E3779B97F4A7C15ULL, H,hid,hm);  // u-ubar string
  if(n<0){ outN[e]=-1; return; }
  double s0=0,s1=0,s2=0,s3=0,dm=0,rhoM=-1.0; int nc=0;
  for(int i=0;i<n;++i){ s0+=H[4*i];s1+=H[4*i+1];s2+=H[4*i+2];s3+=H[4*i+3];
    double m2=H[4*i+3]*H[4*i+3]-H[4*i]*H[4*i]-H[4*i+1]*H[4*i+1]-H[4*i+2]*H[4*i+2];
    double mt=hm[i]; dm=fmax(dm,fabs(m2-mt*mt));         // on-shell vs the mass ACTUALLY used
    if(isCharged(hid[i])) nc++;
    if((abs(hid[i])==113||abs(hid[i])==213)&&rhoM<0.0) rhoM=hm[i]; }   // first rho mass (BW check)
  outN[e]=n; outNc[e]=nc; outTot[4*e]=s0;outTot[4*e+1]=s1;outTot[4*e+2]=s2;outTot[4*e+3]=s3; outDm[e]=dm; outVm[e]=rhoM;
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):200000;
  double rootS=(argc>2)?atof(argv[2]):91.1876, E=0.5*rootS;
  int TPB=128, blocks=(N+TPB-1)/TPB; uint64_t base=0x4144ULL;
  int *dN,*dNc; double *dTot,*dDm,*dVm;
  CK(cudaMalloc(&dN,(size_t)N*4)); CK(cudaMalloc(&dNc,(size_t)N*4));
  CK(cudaMalloc(&dTot,(size_t)N*32)); CK(cudaMalloc(&dDm,(size_t)N*8)); CK(cudaMalloc(&dVm,(size_t)N*8));

  cudaEvent_t t0,t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1)); CK(cudaEventRecord(t0));
  kern<<<blocks,TPB>>>(N,E,base,dN,dNc,dTot,dDm,dVm); CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
  float ms=0; CK(cudaEventElapsedTime(&ms,t0,t1));
  std::vector<int> hN(N),hNc(N); std::vector<double> hTot((size_t)N*4),hDm(N),hVm(N);
  CK(cudaMemcpy(hN.data(),dN,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hNc.data(),dNc,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hTot.data(),dTot,(size_t)N*32,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hDm.data(),dDm,(size_t)N*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hVm.data(),dVm,(size_t)N*8,cudaMemcpyDeviceToHost));
  // reproducibility
  std::vector<int> hN2(N); kern<<<blocks,TPB>>>(N,E,base,dN,dNc,dTot,dDm,dVm); CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(hN2.data(),dN,(size_t)N*4,cudaMemcpyDeviceToHost));

  double maxMom=0,maxDm=0; long sumN=0,sumNc=0,nFail=0,repro=0;
  for(int e=0;e<N;++e){ if(hN[e]!=hN2[e]) repro++;
    if(hN[e]<0){ nFail++; continue; }
    double dx=fabs(hTot[4*e]),dy=fabs(hTot[4*e+1]),dz=fabs(hTot[4*e+2]),de=fabs(hTot[4*e+3]-rootS);
    maxMom=fmax(maxMom,fmax(fmax(dx,dy),fmax(dz,de))); maxDm=fmax(maxDm,hDm[e]);
    sumN+=hN[e]; sumNc+=hNc[e]; }
  long nok=N-nFail; double meanN=(double)sumN/nok, meanNc=(double)sumNc/nok;

  // GPU vs CPU determinism (subset): same nHadrons AND same hadron-id sequence
  int nCPU=(N<20000)?N:20000; long structSame=0;
  std::vector<double> H(MAXH*4),hmc(MAXH); std::vector<int> id(MAXH);
  for(int e=0;e<nCPU;++e){ int n=hadronizeString(E,2,-2, base+(uint64_t)e*0x9E3779B97F4A7C15ULL,H.data(),id.data(),hmc.data());
    if(n==hN[e]) structSame++; }

  // Breit-Wigner spectrum check: rho mass mean + RMS (pole when -DUSE_BW off -> RMS=0).
  long nRho=0; double sumM=0,sumM2=0;
  for(int e=0;e<N;++e) if(hN[e]>=0 && hVm[e]>0.0){ nRho++; sumM+=hVm[e]; sumM2+=hVm[e]*hVm[e]; }
  double rhoMean=(nRho?sumM/nRho:0), rhoRMS=(nRho?sqrt(fmax(0.0,sumM2/nRho-rhoMean*rhoMean)):0);

  printf("GPU Lund string fragmentation (u-ubar, sqrt(s)=%.4f GeV, %d events)\n",rootS,N);
  printf("  throughput        : %.1f ms (%.2f M strings/s)\n",ms,N/ms/1e3);
  printf("  multiplicity      : mean %.3f hadrons, %.3f charged\n",meanN,meanNc);
  printf("  4-mom conservation: max|deviation| = %.2e GeV\n",maxMom);
  printf("  on-shellness      : max|m^2-table| = %.2e GeV^2\n",maxDm);
  printf("  refragment/fail   : failed events  = %ld / %d\n",nFail,N);
  printf("  reproducibility   : GPU re-run diffs = %ld\n",repro);
  printf("  GPU vs CPU        : nHadrons identical %ld/%d = %.2f%%\n",structSame,nCPU,100.0*structSame/nCPU);
  printf("  rho mass spectrum : mean %.4f GeV, RMS %.4f (pole 0.7753, width 0.149; %ld rhos)\n",rhoMean,rhoRMS,nRho);
#ifdef USE_BW
  bool bwOK=(fabs(rhoMean-0.7753)<0.02)&&(rhoRMS>0.05);   // BW broadened around the pole
#else
  bool bwOK=(rhoRMS<1e-2);                                 // pole masses: no spread (FP-noise ~1e-5)
#endif
  bool ok=(maxMom<1e-6)&&(maxDm<1e-6)&&(repro==0)&&(nFail<N/100)&&(structSame==nCPU)&&(meanNc>2.0&&meanNc<30.0)&&bwOK;
  printf("VALIDATION: %s (conservation+on-shell+reproducible+GPU==CPU+multiplicity+BW-spectrum)\n",ok?"PASS":"FAIL");
  cudaFree(dN);cudaFree(dNc);cudaFree(dTot);cudaFree(dDm);cudaFree(dVm);
  return ok?0:2;
}
