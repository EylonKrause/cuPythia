// cuPythia pipeline stage 3 — PHYSICAL final-state (timelike) dipole shower on the GPU.
//
// This is the headline novel piece: a Sudakov-veto parton shower running one event per
// GPU thread (the GAPS decomposition, arXiv:2403.08692 / 2511.19633, Seymour & Sule),
// with the splitting kernels, running-alpha_s trial generation, z-sampling and exact
// local-dipole RECOIL kinematics ported from Pythia 8.317 SimpleTimeShower.
//
// Physics scope (honest): e+e- -> Z -> q qbar at the Z pole, FINAL-STATE radiation only,
// massless partons, splittings q->qg and g->gg (g->qqbar omitted — flagged below), 1-loop
// running alpha_s with FIXED n_f=5 (no flavour-threshold matching — a deliberate, labelled
// simplification). The colour chain is large-N_c planar: partons are kept in colour order
// (q ... gluons ... qbar) so a dipole is simply an adjacent pair and an emission is an
// insertion — exactly the dipole-shower picture GAPS uses.
//
// Validation: (1) GPU run twice -> bit-identical (counter-RNG reproducibility);
//             (2) exact 4-momentum conservation + on-shellness of every final parton;
//             (3) GPU vs an IDENTICAL CPU port -> same mean multiplicity + per-event
//                 bit-identity fraction (FP transcendental ULPs can flip a veto decision).
// Rivet observables (thrust, Durham jet rates) vs Pythia are the next validation layer.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o shower_fsr shower_fsr.cu
// Run:   ./shower_fsr [nEvents=200000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

static const double MZ     = 91.1876;       // Z mass / CM energy
static const double EBEAM  = 0.5*MZ;        // each initial quark energy
static const double PT2MIN = 0.5*0.5;       // (TimeShower:pTmin = 0.5 GeV)^2 cutoff
static const double ASMZ   = 0.1365;        // alpha_s(M_Z) (TimeShower:alphaSvalue)
#define MAXP 64

// 1-loop running alpha_s referenced to alpha_s(M_Z)=0.1365, with FLAVOUR-THRESHOLD
// matching (n_f = 5,4,3 across m_b, m_c) — closes one residual vs Pythia, which fixed
// n_f=5 did not. 1/alpha runs continuously; below each threshold beta0 grows so alpha
// rises faster (the physically larger low-scale coupling Pythia also uses).
__host__ __device__ inline double alphaS(double mu2){
  const double mc2=1.5*1.5, mb2=4.8*4.8, mZ2=MZ*MZ;
  const double b5=(33.0-2.0*5.0)/(12.0*M_PI), b4=(33.0-2.0*4.0)/(12.0*M_PI), b3=(33.0-2.0*3.0)/(12.0*M_PI);
  double inv=1.0/ASMZ;                         // 1/alpha at M_Z (n_f=5)
  if(mu2>=mb2){ inv+=b5*log(mu2/mZ2); }
  else { inv+=b5*log(mb2/mZ2);                 // n_f=5: M_Z -> m_b
    if(mu2>=mc2){ inv+=b4*log(mu2/mb2); }      // n_f=4: m_b -> mu
    else { inv+=b4*log(mc2/mb2)+b3*log(mu2/mc2); } } // n_f=4: m_b->m_c, n_f=3: m_c->mu
  double a=1.0/inv;
  return (a>0.0 && a<10.0)? a : 10.0;          // guard the Landau region (mu2 stays >= PT2MIN)
}
__host__ __device__ inline double mass2(const double* a,const double* b){
  double e=a[3]+b[3],x=a[0]+b[0],y=a[1]+b[1],z=a[2]+b[2]; return e*e-x*x-y*y-z*z;
}
__host__ __device__ inline void cp4(double* d,const double* s){ d[0]=s[0];d[1]=s[1];d[2]=s[2];d[3]=s[3]; }

// Active Lorentz boost of q by velocity beta=(ex,ey,ez), Lorentz factor gamma.
__host__ __device__ inline void boostBy(const double* q,double ex,double ey,double ez,double gamma,double* o){
  double bdq=ex*q[0]+ey*q[1]+ez*q[2], e2=ex*ex+ey*ey+ez*ez;
  double k=(e2>1e-18)?((gamma-1.0)*bdq/e2+gamma*q[3]):0.0;
  o[0]=q[0]+k*ex; o[1]=q[1]+k*ey; o[2]=q[2]+k*ez; o[3]=gamma*(q[3]+bdq);
}
// Rotate the spatial part of v so that the +z axis maps to unit vector n.
__host__ __device__ inline void rotZto(double nx,double ny,double nz,double* v){
  double cz=(nz>1?1:(nz<-1?-1:nz)); double th=acos(cz), ph=atan2(ny,nx);
  double ct=cos(th),st=sin(th),cp=cos(ph),sp=sin(ph);
  double x=v[0],y=v[1],z=v[2];
  double x1=x*ct+z*st, y1=y, z1=-x*st+z*ct;
  v[0]=x1*cp-y1*sp; v[1]=x1*sp+y1*cp; v[2]=z1;
}

// Pythia SimpleTimeShower::branch kinematics (massless rad/rec/emt): build the post-emission
// radiator/emitted/recoiler momenta in the lab from (pT2,z,phi). Returns false if unphysical.
__host__ __device__ inline bool doKin(const double* R,const double* C,double pT2,double z,double phi,
                                      double* oR,double* oE,double* oC){
  double m2Dip=mass2(R,C); if(m2Dip<=0) return false; double mDip=sqrt(m2Dip);
  double m2=pT2/(z*(1.0-z));               if(m2>=m2Dip) return false;
  double eRPE =0.5*(m2Dip+m2)/mDip, pzRPE=0.5*(m2Dip-m2)/mDip;
  double pT2c =m2*(eRPE*eRPE*z*(1.0-z)-0.25*m2)/(pzRPE*pzRPE); if(pT2c<0) return false;
  double pTc=sqrt(pT2c);
  double pzRad=(eRPE*eRPE*z-0.5*m2)/pzRPE, pzEmt=(eRPE*eRPE*(1.0-z)-0.5*m2)/pzRPE;
  double cph=cos(phi),sph=sin(phi);
  double qR[4]={ pTc*cph, pTc*sph, pzRad, sqrt(pTc*pTc+pzRad*pzRad)};
  double qE[4]={-pTc*cph,-pTc*sph, pzEmt, sqrt(pTc*pTc+pzEmt*pzEmt)};
  double qC[4]={0.0,0.0,-pzRPE, pzRPE};
  double PE=R[3]+C[3], gamma=PE/mDip, bx=(R[0]+C[0])/PE,by=(R[1]+C[1])/PE,bz=(R[2]+C[2])/PE;
  double Rcm[4]; boostBy(R,-bx,-by,-bz,gamma,Rcm);   // radiator direction in the dipole CM
  double nn=sqrt(Rcm[0]*Rcm[0]+Rcm[1]*Rcm[1]+Rcm[2]*Rcm[2]);
  double nx,ny,nz; if(nn<1e-12){nx=0;ny=0;nz=1;}else{nx=Rcm[0]/nn;ny=Rcm[1]/nn;nz=Rcm[2]/nn;}
  rotZto(nx,ny,nz,qR); rotZto(nx,ny,nz,qE); rotZto(nx,ny,nz,qC);
  boostBy(qR,bx,by,bz,gamma,oR); boostBy(qE,bx,by,bz,gamma,oE); boostBy(qC,bx,by,bz,gamma,oC);
  return true;
}

// Shower one event into the local arrays P (px,py,pz,E per parton) and id; return nPartons.
__host__ __device__ inline int showerEvent(double* P,int* id,uint64_t ctr){
  P[0]=0;P[1]=0;P[2]= EBEAM;P[3]=EBEAM; id[0]= 1;   // q   (colour end of chain)
  P[4]=0;P[5]=0;P[6]=-EBEAM;P[7]=EBEAM; id[1]=-1;   // qbar(anticolour end)
  int n=2; double pT2=0.25*MZ*MZ;                    // pT2begDip = 0.25 m2Dip
  for(int step=0; step<MAXP; ++step){
    double bestT=PT2MIN, bestZ=0; int bRad=-1,bRec=-1;
    // enumerate colour (i->i+1) and anticolour (i->i-1) dipole ends
    for(int i=0;i<n;++i){
      for(int side=0;side<2;++side){
        int rec; bool valid;
        if(side==0){ rec=i+1; valid=(rec<n)&&(id[i]>0||id[i]==21); }   // colour end
        else       { rec=i-1; valid=(rec>=0)&&(id[i]<0||id[i]==21); }  // anticolour end
        if(!valid) continue;
        double m2Dip=mass2(P+4*i,P+4*rec); if(m2Dip<=4.0*PT2MIN) continue;
        bool isQ=(id[i]!=21); double colFac=isQ?(4.0/3.0):(3.0/2.0);
        double zc=0.5-sqrt(fmax(0.0,0.25-PT2MIN/m2Dip)); if(zc<1e-10) zc=PT2MIN/m2Dip;
        if(zc>0.499) continue; double zmaxc=1.0-zc;
        double amax=alphaS(PT2MIN);
        double Iover=(amax/(2.0*M_PI))*colFac*2.0*log(zmaxc/zc); if(Iover<=0) continue;
        double t=pT2;
        for(int it=0; it<20000; ++it){
          double R1=u01(splitmix64(ctr++)); t*=pow(R1,1.0/Iover);
          if(t<PT2MIN) break;
          double R2=u01(splitmix64(ctr++));
          double zz=1.0-(1.0-zc)*pow(zc/(1.0-zc),R2);
          double R3=u01(splitmix64(ctr++));
          double zmn=0.5-sqrt(fmax(0.0,0.25-t/m2Dip)), zmx=1.0-zmn;
          double m2v=t/(zz*(1.0-zz));
          if(zz>zmn && zz<zmx && m2v<m2Dip){
            double dal=zz*(1.0-zz)*(m2Dip+m2v)*(m2Dip+m2v);          // Dalitz bound (m2Rec=0)
            if(m2v*m2Dip<dal){
              double w=isQ?0.5*(1.0+zz*zz):0.5*(1.0+zz*zz*zz);       // (1+z^2)/2 or (1+z^3)/2
              if(R3<(alphaS(t)/amax)*w){
                if(t>bestT){ bestT=t; bestZ=zz; bRad=i; bRec=rec; }
                break;
              }
            }
          }
          // else veto: t already lowered, continue
        }
      }
    }
    if(bRad<0) break;                  // no end emitted above the cutoff -> shower done
    pT2=bestT;
    double phi=2.0*M_PI*u01(splitmix64(ctr++));
    double oR[4],oE[4],oC[4];
    if(!doKin(P+4*bRad,P+4*bRec,bestT,bestZ,phi,oR,oE,oC)) continue;
    if(n>=MAXP) break;
    if(bRec==bRad+1){                  // colour end: gluon between rad and rad+1
      for(int k=n;k>bRad+1;--k){ cp4(P+4*k,P+4*(k-1)); id[k]=id[k-1]; }
      cp4(P+4*bRad,oR);
      cp4(P+4*(bRad+1),oE); id[bRad+1]=21;
      cp4(P+4*(bRad+2),oC);
    } else {                           // anticolour end: gluon between rad-1 and rad
      for(int k=n;k>bRad;--k){ cp4(P+4*k,P+4*(k-1)); id[k]=id[k-1]; }
      cp4(P+4*(bRad-1),oC);
      cp4(P+4*bRad,oE); id[bRad]=21;
      cp4(P+4*(bRad+1),oR);
    }
    n++;
  }
  return n;
}

// Event thrust T = max_n  sum_i |p_i . n| / sum_i |p_i|, by the standard iterative
// fixed-point (n -> sum_i sign(p_i.n) p_i, renormalise) with a few seed axes.
__host__ __device__ inline double thrust(const double* P,int n){
  double psum=0; for(int i=0;i<n;++i) psum+=sqrt(P[4*i]*P[4*i]+P[4*i+1]*P[4*i+1]+P[4*i+2]*P[4*i+2]);
  if(psum<=0) return 1.0;
  double Tbest=0;
  for(int s=0;s<n;++s){           // seed from every particle direction (robust for multi-jet events)
    double ax=P[4*s],ay=P[4*s+1],az=P[4*s+2]; double a=sqrt(ax*ax+ay*ay+az*az);
    if(a<1e-12) continue; ax/=a;ay/=a;az/=a;
    for(int it=0;it<20;++it){ double nx=0,ny=0,nz=0;
      for(int i=0;i<n;++i){ double d=P[4*i]*ax+P[4*i+1]*ay+P[4*i+2]*az; double sg=(d>=0)?1.0:-1.0;
        nx+=sg*P[4*i];ny+=sg*P[4*i+1];nz+=sg*P[4*i+2]; }
      double nn=sqrt(nx*nx+ny*ny+nz*nz); if(nn<1e-12) break; ax=nx/nn;ay=ny/nn;az=nz/nn; }
    double num=0; for(int i=0;i<n;++i) num+=fabs(P[4*i]*ax+P[4*i+1]*ay+P[4*i+2]*az);
    Tbest=fmax(Tbest,num/psum);
  }
  return Tbest;
}

__global__ void showerKernel(int nEvt,uint64_t base,int* outN,double* outTot,double* outM2,double* outThr){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=nEvt) return;
  double P[MAXP*4]; int id[MAXP];
  int n=showerEvent(P,id, base + (uint64_t)e*0x9E3779B97F4A7C15ULL);
  double s0=0,s1=0,s2=0,s3=0,mm=0;
  for(int i=0;i<n;++i){ s0+=P[4*i];s1+=P[4*i+1];s2+=P[4*i+2];s3+=P[4*i+3];
    double m2=P[4*i+3]*P[4*i+3]-P[4*i]*P[4*i]-P[4*i+1]*P[4*i+1]-P[4*i+2]*P[4*i+2];
    mm=fmax(mm,fabs(m2)); }
  outN[e]=n; outTot[4*e]=s0;outTot[4*e+1]=s1;outTot[4*e+2]=s2;outTot[4*e+3]=s3; outM2[e]=mm;
  outThr[e]=thrust(P,n);
}

int main(int argc,char**argv){
  int nEvt=(argc>1)?atoi(argv[1]):200000;
  int TPB=128, blocks=(nEvt+TPB-1)/TPB;        // GAPS-optimal 128 threads/block
  uint64_t base=0x5110UL;

  int *dN; double *dTot,*dM2,*dThr;
  CK(cudaMalloc(&dN,(size_t)nEvt*4)); CK(cudaMalloc(&dTot,(size_t)nEvt*32));
  CK(cudaMalloc(&dM2,(size_t)nEvt*8)); CK(cudaMalloc(&dThr,(size_t)nEvt*8));

  cudaEvent_t t0,t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
  CK(cudaEventRecord(t0));
  showerKernel<<<blocks,TPB>>>(nEvt,base,dN,dTot,dM2,dThr);
  CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
  float ms=0; CK(cudaEventElapsedTime(&ms,t0,t1));

  std::vector<int> hN(nEvt); std::vector<double> hTot((size_t)nEvt*4),hM2(nEvt),hThr(nEvt);
  CK(cudaMemcpy(hN.data(),dN,(size_t)nEvt*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hTot.data(),dTot,(size_t)nEvt*32,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hM2.data(),dM2,(size_t)nEvt*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hThr.data(),dThr,(size_t)nEvt*8,cudaMemcpyDeviceToHost));

  // (1) reproducibility: a second identical launch must be bit-identical.
  std::vector<int> hN2(nEvt); std::vector<double> hTot2((size_t)nEvt*4);
  showerKernel<<<blocks,TPB>>>(nEvt,base,dN,dTot,dM2,dThr); CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(hN2.data(),dN,(size_t)nEvt*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hTot2.data(),dTot,(size_t)nEvt*32,cudaMemcpyDeviceToHost));
  long reproDiff=0; for(int e=0;e<nEvt;++e){ if(hN[e]!=hN2[e]) reproDiff++;
    for(int k=0;k<4;++k) if(hTot[4*e+k]!=hTot2[4*e+k]) {reproDiff++;break;} }

  // (2) physics: 4-momentum conservation + on-shellness; multiplicity stats.
  double maxMomViol=0, maxM2=0; long sumN=0; int minN=1<<30,maxN=0;
  for(int e=0;e<nEvt;++e){ double dx=fabs(hTot[4*e]),dy=fabs(hTot[4*e+1]),dz=fabs(hTot[4*e+2]),de=fabs(hTot[4*e+3]-MZ);
    double v=fmax(fmax(dx,dy),fmax(dz,de)); maxMomViol=fmax(maxMomViol,v);
    maxM2=fmax(maxM2,hM2[e]); sumN+=hN[e]; if(hN[e]<minN)minN=hN[e]; if(hN[e]>maxN)maxN=hN[e]; }
  double meanN=(double)sumN/nEvt;

  // (3) identical CPU port over a subset. The shower CONTROL FLOW (multiplicity, i.e. every
  //     accept/veto decision) must be bit-identical to the GPU; the per-event summed momenta
  //     then agree only to GPU-vs-CPU IEEE transcendental accumulation (never bit-identical).
  int nCPU=(nEvt<20000)?nEvt:20000; long sumNcpu=0,structSame=0; double maxMomRel=0;
  std::vector<double> P(MAXP*4); std::vector<int> id(MAXP);
  for(int e=0;e<nCPU;++e){ int n=showerEvent(P.data(),id.data(), base+(uint64_t)e*0x9E3779B97F4A7C15ULL);
    sumNcpu+=n; double s0=0,s1=0,s2=0,s3=0;
    for(int i=0;i<n;++i){ s0+=P[4*i];s1+=P[4*i+1];s2+=P[4*i+2];s3+=P[4*i+3]; }
    if(n==hN[e]){ structSame++;
      double d=fmax(fmax(fabs(s0-hTot[4*e]),fabs(s1-hTot[4*e+1])),fmax(fabs(s2-hTot[4*e+2]),fabs(s3-hTot[4*e+3])));
      maxMomRel=fmax(maxMomRel,d/MZ); } }
  double meanNcpu=(double)sumNcpu/nCPU;

  // thrust observable: mean(1-T) and a normalised (1-T) histogram dumped for Pythia comparison.
  const int NB=20; const double TMAX=0.5; long hist[NB]={0}; double sum1mT=0;
  for(int e=0;e<nEvt;++e){ double omt=1.0-hThr[e]; sum1mT+=omt;
    int b=(int)(omt/TMAX*NB); if(b<0)b=0; if(b>=NB)b=NB-1; hist[b]++; }
  double mean1mT=sum1mT/nEvt;
  FILE* fh=fopen("thrust_gpu.dat","w");
  if(fh){ fprintf(fh,"# (1-T)_low  (1-T)_high  normalised_density   [cuPythia GPU FSR shower, %d evts]\n",nEvt);
    for(int b=0;b<NB;++b) fprintf(fh,"%.4f %.4f %.6e\n",b*TMAX/NB,(b+1)*TMAX/NB,hist[b]/((double)nEvt*(TMAX/NB)));
    fclose(fh); }

  printf("FSR dipole shower on GPU (e+e- -> Z -> q qbar, sqrt(s)=%.4f GeV, %d events)\n",MZ,nEvt);
  printf("  throughput        : %.2f ms  (%.2f M events/s)\n", ms, nEvt/ms/1e3);
  printf("  multiplicity      : mean %.3f partons  (min %d, max %d)\n", meanN, minN, maxN);
  printf("  4-mom conservation: max|deviation| = %.2e GeV\n", maxMomViol);
  printf("  on-shellness      : max|p^2|        = %.2e GeV^2\n", maxM2);
  printf("  thrust            : <1-T> = %.4f  (histogram -> thrust_gpu.dat)\n", mean1mT);
  printf("  reproducibility   : GPU re-run diffs = %ld  (counter-RNG)\n", reproDiff);
  printf("  GPU vs CPU port   : control-flow bit-identical %ld/%d = %.2f%%  (mean mult %.3f vs %.3f)\n",
         structSame, nCPU, 100.0*structSame/nCPU, meanN, meanNcpu);
  printf("                      momenta agree to %.2e (GPU/CPU IEEE transcendental accumulation)\n", maxMomRel);
  bool ok = (maxMomViol<1e-5) && (maxM2<1e-3) && (reproDiff==0) &&
            (meanN>2.5&&meanN<25.0) && (structSame==nCPU) && (maxMomRel<1e-6) &&
            (mean1mT>0.01 && mean1mT<0.30);
  printf("VALIDATION: %s (momentum+on-shell+reproducible+CPU-agreement+thrust-sane)\n", ok?"PASS":"FAIL");
  cudaFree(dN);cudaFree(dTot);cudaFree(dM2);cudaFree(dThr);
  return ok?0:2;
}
