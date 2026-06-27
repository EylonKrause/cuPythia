// Standalone validation of the Lund baryon FLAVOUR logic (pickFlavBreakMR / combineBaryonMR /
// combineHadronMR), decoupled from the heavy hadronize_mr kernel so it compiles fast. Simulates the
// flavour side of string fragmentation (a q ... qbar string broken from both ends + finalTwo) and
// checks: (1) electric-charge conservation, (2) baryon-number conservation -- both must be EXACT 0
// per event (an e+e- -> q qbar string is charge- and baryon-number-neutral); plus the baryon rates.
// The functions below are copied verbatim from hadronize_mr.cu (the physics under test).
//
// Build: nvcc -O2 -std=c++17 -arch=sm_120 -o baryon_test baryon_test.cu
// Run:   ./baryon_test [nEvents=2000000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"
#define CK(c) do{cudaError_t e_=(c); if(e_!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e_),__LINE__);return 1;}}while(0)

// ---- flavour layer copied from hadronize_mr.cu (params + meson combine + baryon mechanism) ----
static const double H_PROBSTOUD=0.217, H_MUDV=0.50, H_SV=0.55, H_THPS=-15.0, H_THV=36.0;
static const double H_ETASUP=0.60, H_ETAPSUP=0.12;
__host__ __device__ inline int pickFlavMR(uint64_t& c){ double r=u01(splitmix64(c++))*(2.0+H_PROBSTOUD); return (r<1.0)?1:((r<2.0)?2:3); }
__host__ __device__ inline int combineMesonMR(int id1,int id2,uint64_t& c){
  int a1=abs(id1),a2=abs(id2),idMax=(a1>a2)?a1:a2,idMin=(a1<a2)?a1:a2;
  int flav=(idMax<3)?0:idMax-2; double mv=(flav==0)?H_MUDV:H_SV;
  double rs=(1.0+mv)*u01(splitmix64(c++)); int spin=0; rs-=1.0; if(rs>0.0) spin=1;
  int code=(spin==0)?1:3, idM=100*idMax+10*idMin+code;
  if(idMax!=idMin){ int sg=(idMax%2==0)?1:-1; if((idMax==a1&&id1<0)||(idMax==a2&&id2<0)) sg=-sg; idM*=sg; }
  else { double al=((spin==0)?(90.0-(H_THPS+54.7)):(H_THV+54.7))*M_PI/180.0;
    double m1,m2; if(flav==0){m1=0.5;m2=0.5*(1.0+sin(al)*sin(al));}else{m1=0.0;m2=cos(al)*cos(al);}
    double rm=u01(splitmix64(c++)); if(rm<m1)idM=110;else if(rm<m2)idM=220;else idM=330; idM+=code;
    if(idM==221&&H_ETASUP<u01(splitmix64(c++)))return 0; if(idM==331&&H_ETAPSUP<u01(splitmix64(c++)))return 0; }
  return idM;
}
static const double H_PROBQQTOQ=0.081, H_PROBSQTOQQ=0.915, H_PROBQQ1TOQQ0=0.0275, H_DECUPLETSUP=1.0;
__host__ __device__ inline double bCGoct(int i){ const double v[6]={0.75,0.5,0.0,0.1667,0.0833,0.1667}; return v[i]; }
__host__ __device__ inline double bCGdec(int i){ const double v[6]={0.0,0.0,1.0,0.3333,0.6667,0.3333}; return v[i]; }
__host__ __device__ inline double bCGsum(int i){ return bCGoct(i)+H_DECUPLETSUP*bCGdec(i); }
__host__ __device__ inline double bCGmax(int i){
  double m01=fmax(bCGsum(0),bCGsum(1)), m23=fmax(bCGsum(2),bCGsum(3)), m45=fmax(bCGsum(4),bCGsum(5));
  return (i<2)?m01:((i<4)?m23:m45); }
__host__ __device__ inline int pickDiqFlav(uint64_t& c){
  double r=(2.0+H_PROBSQTOQQ*H_PROBSTOUD)*u01(splitmix64(c++)); return (r<1.0)?1:((r<2.0)?2:3); }
__host__ __device__ inline int pickFlavBreakMR(int oldFlav,uint64_t& c){
  int idOld=abs(oldFlav); bool doOldB=(idOld>1000);
  bool roll=((1.0+H_PROBQQTOQ)*u01(splitmix64(c++))>1.0);
  bool doNewB=(!doOldB)&&(idOld<4)&&roll;
  if(!doNewB){ int q=pickFlavMR(c); if((oldFlav>0&&oldFlav<9)||oldFlav<-1000) q=-q; return q; }
  int idP=pickDiqFlav(c), idV=pickDiqFlav(c);
  if(idP<3&&idV<3){ idV=idP; if(u01(splitmix64(c++))>0.5) idV=3-idP; } else (void)u01(splitmix64(c++));
  int spin=3; if(idV!=idP){ if((1.0+3.0*H_PROBQQ1TOQQ0)*u01(splitmix64(c++))<1.0) spin=1; } else (void)u01(splitmix64(c++));
  int dq=1000*((idV>idP)?idV:idP)+100*((idV<idP)?idV:idP)+spin;
  if((oldFlav<0&&oldFlav>-9)||oldFlav>1000) dq=-dq; return dq; }
__host__ __device__ inline int combineBaryonMR(int id1,int id2,uint64_t& c){
  int a1=abs(id1),a2=abs(id2),idMax=(a1>a2)?a1:a2,idMin=(a1<a2)?a1:a2;
  int idQQ1=idMax/1000, idQQ2=(idMax/100)%10, spinQQ=idMax%10;
  int sf=spinQQ-1; if(sf==2&&idQQ1!=idQQ2) sf=4; if(idMin!=idQQ1&&idMin!=idQQ2) sf++;
  if(sf<0||sf>5) return 0;
  if(bCGsum(sf) < u01(splitmix64(c++))*bCGmax(sf)) return 0;
  int o1=(idMin>idQQ1)?((idMin>idQQ2)?idMin:idQQ2):((idQQ1>idQQ2)?idQQ1:idQQ2);
  int o3=(idMin<idQQ1)?((idMin<idQQ2)?idMin:idQQ2):((idQQ1<idQQ2)?idQQ1:idQQ2);
  int o2=idMin+idQQ1+idQQ2-o1-o3;
  int spinBar=(bCGsum(sf)*u01(splitmix64(c++)) < bCGoct(sf))?2:4;
  bool lam=false;
  if(spinBar==2&&o1>o2&&o2>o3){ lam=(spinQQ==1);
    if(o1!=idMin&&spinQQ==1) lam=(u01(splitmix64(c++))<0.25);
    else if(o1!=idMin)       lam=(u01(splitmix64(c++))<0.75);
    else (void)u01(splitmix64(c++)); }
  else (void)u01(splitmix64(c++));
  int idB=lam?(1000*o1+100*o3+10*o2+spinBar):(1000*o1+100*o2+10*o3+spinBar);
  return (id1>0)?idB:-idB; }
__host__ __device__ inline int combineHadronMR(int id1,int id2,uint64_t& c){
  if(abs(id1)>1000||abs(id2)>1000) return combineBaryonMR(id1,id2,c);
  return combineMesonMR(id1,id2,c); }

// 3*charge of a hadron (from PDG digits); meson q1 q2bar; baryon 3 quarks. Used only by the checker.
__host__ __device__ inline int chq3(int q){ q=abs(q)%10; return (q==2||q==4)?2:-1; }   // u,c=+2/3; d,s,b=-1/3
__host__ __device__ inline int charge3(int pid){   // 3*electric charge
  int a=abs(pid), s=(pid>0)?1:-1;
  if(a>1000){ int q1=(a/1000),q2=(a/100)%10,q3=(a/10)%10; return s*(chq3(q1)+chq3(q2)+chq3(q3)); } // baryon: sum of 3 quarks
  // charged mesons: the +PDG code is the +1 particle (pi+,rho+,K+,K*+,D+,D*+,Ds+,Ds*+,B+,B*+,Bc+,Bc*+);
  // all diagonal/neutral mesons (pi0,K0,eta,eta',rho0,omega,phi,K*0,K_L,K_S,B0,Bs0,...) and gamma = 0.
  if(a==211||a==213||a==321||a==323||a==411||a==413||a==431||a==433||a==521||a==523||a==541||a==543) return 3*s;
  return 0;
}
__host__ __device__ inline int baryonNum(int pid){ return (abs(pid)>1000)?((pid>0)?1:-1):0; }

// Simulate the FLAVOUR side of fragmenting one q...qbar string: alternate breaks off the two ends,
// then a finalTwo, tallying charge3 and baryon number. Returns packed (sumQ3, sumB, nBaryon, codeOfP/L...).
__host__ __device__ inline void fragFlav(uint64_t ctr,int& sumQ3,int& sumB,int& nBar,int& nP,int& nLam,int& nMesonBad){
  int posF=1, negF=-1; sumQ3=0; sumB=0; nBar=0; nP=0; nLam=0; nMesonBad=0;
  int nBreak=8+(int)(8.0*u01(splitmix64(ctr++)));   // 8..15 breaks then finalTwo (a typical string)
  for(int s=0;s<nBreak;++s){
    bool fromPos=(u01(splitmix64(ctr++))<0.5);
    int e = fromPos?posF:negF;
    int pdg=0, fSel=0;
    for(int tr=0; tr<20 && pdg==0; ++tr){ fSel=pickFlavBreakMR(e,ctr); pdg=combineHadronMR(e,fSel,ctr); }
    if(pdg==0) continue;
    sumQ3+=charge3(pdg); sumB+=baryonNum(pdg); if(abs(pdg)>1000)nBar++;
    if(abs(pdg)==2212)nP++; if(abs(pdg)==3122)nLam++;
    if(charge3(pdg)>10||charge3(pdg)<-10) nMesonBad++;       // sanity: |charge|<=2
    int ne=-fSel; if(fromPos)posF=ne; else negF=ne;
  }
  // finalTwo: single q-qbar pair joins the two ends.
  int p1=0,p2=0;
  for(int tr=0; tr<40 && (p1==0||p2==0); ++tr){ int fN=pickFlavMR(ctr);
    p1=combineHadronMR(posF,-fN,ctr); p2=combineHadronMR(fN,negF,ctr); }
  if(p1){ sumQ3+=charge3(p1); sumB+=baryonNum(p1); if(abs(p1)>1000)nBar++; if(abs(p1)==2212)nP++; if(abs(p1)==3122)nLam++; }
  if(p2){ sumQ3+=charge3(p2); sumB+=baryonNum(p2); if(abs(p2)>1000)nBar++; if(abs(p2)==2212)nP++; if(abs(p2)==3122)nLam++; }
}

__global__ void kern(int N,uint64_t base,int* oQ,int* oB,int* oBar,int* oP,int* oL){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  int q,b,nb,np,nl,bad; fragFlav(base+(uint64_t)e*0x9E3779B97F4A7C15ULL,q,b,nb,np,nl,bad);
  oQ[e]=q; oB[e]=b; oBar[e]=nb; oP[e]=np; oL[e]=nl;
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):2000000; int TPB=128, blocks=(N+TPB-1)/TPB; uint64_t base=0xBA47;
  int *dQ,*dB,*dBar,*dP,*dL;
  CK(cudaMalloc(&dQ,(size_t)N*4));CK(cudaMalloc(&dB,(size_t)N*4));CK(cudaMalloc(&dBar,(size_t)N*4));CK(cudaMalloc(&dP,(size_t)N*4));CK(cudaMalloc(&dL,(size_t)N*4));
  kern<<<blocks,TPB>>>(N,base,dQ,dB,dBar,dP,dL); CK(cudaDeviceSynchronize());
  std::vector<int> hQ(N),hB(N),hBar(N),hP(N),hL(N);
  CK(cudaMemcpy(hQ.data(),dQ,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hB.data(),dB,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hBar.data(),dBar,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hP.data(),dP,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hL.data(),dL,(size_t)N*4,cudaMemcpyDeviceToHost));
  long badQ=0,badB=0,sumBar=0,sumP=0,sumL=0;
  for(int e=0;e<N;++e){ if(hQ[e]!=0)badQ++; if(hB[e]!=0)badB++; sumBar+=hBar[e]; sumP+=hP[e]; sumL+=hL[e]; }
  printf("Lund baryon flavour test (%d strings, q..qbar broken both ends + finalTwo)\n",N);
  printf("  charge non-conserving strings        : %ld / %d\n",badQ,N);
  printf("  baryon-number non-conserving strings : %ld / %d\n",badB,N);
  printf("  baryons/string %.3f   p+pbar/string %.3f   Lambda+Lbar/string %.3f\n",(double)sumBar/N,(double)sumP/N,(double)sumL/N);
  bool ok=(badQ==0 && badB==0);
  printf("BARYON-FLAVOUR VALIDATION: %s (charge + baryon number conserved on every string)\n",ok?"PASS":"FAIL");
  cudaFree(dQ);cudaFree(dB);cudaFree(dBar);cudaFree(dP);cudaFree(dL);
  return ok?0:2;
}
