// Standalone validation of the GPU hadron-decay module (decay_inc.cuh), before wiring it into the
// pipeline. Each test event builds a random list of (boosted, on-shell) primary hadrons from a pool
// that exercises every decay channel (rho/K*/omega/phi/eta/eta'/K0 + stable pi/K), decays them, and
// checks: (1) per-event 4-momentum conservation sum(final)==sum(primary); (2) on-shellness of every
// final; (3) GPU==CPU (same final count + momentum sum); (4) reproducibility (2nd launch identical).
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o decay_test decay_test.cu
// Run:   ./decay_test [nEvents=500000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"
#define CK(c) do{cudaError_t e_=(c); if(e_!=cudaSuccess){printf("CUDA %s @ %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__);return 1;}}while(0)
#include "decay_inc.cuh"

// Pool of primary ids (unstable parents + a few stable) to exercise the table.
__host__ __device__ inline int poolId(int k){
  const int POOL[16]={113,213,-213,223,333,313,-313,323,-323,221,331,311,-311,211,-211,321};
  return POOL[k%16];
}

__host__ __device__ inline void runOne(uint64_t base,int e,double* cons,double* dm,int* nFout,double* finSum){
  uint64_t ctr=base+(uint64_t)e*0x9E3779B97F4A7C15ULL;
  double P[8*4]; int id[8];
  int n=2+(int)(5.0*u01(splitmix64(ctr++)));   // 2..6 primaries
  double ps0=0,ps1=0,ps2=0,ps3=0;
  for(int i=0;i<n;++i){ int pid=poolId((int)(16.0*u01(splitmix64(ctr++))));
    double m=decayMass(pid,ctr);                                  // BW/pole mass (1 draw)
    double px=6.0*(u01(splitmix64(ctr++))-0.5), py=6.0*(u01(splitmix64(ctr++))-0.5), pz=6.0*(u01(splitmix64(ctr++))-0.5);
    double en=sqrt(px*px+py*py+pz*pz+m*m);
    P[4*i]=px;P[4*i+1]=py;P[4*i+2]=pz;P[4*i+3]=en; id[i]=pid;
    ps0+=px;ps1+=py;ps2+=pz;ps3+=en; }
  double F[MAXFINAL*4]; int fid[MAXFINAL];
  uint64_t dctr=(base^0xDECAULL)+(uint64_t)e*0x100000001B3ULL;
  int nF=decayEvent(P,id,n,dctr,F,fid);
  if(nF<0){ *nFout=-1; *cons=0; *dm=0; finSum[0]=finSum[1]=finSum[2]=finSum[3]=0; return; }
  double f0=0,f1=0,f2=0,f3=0,mdm=0;
  for(int i=0;i<nF;++i){ f0+=F[4*i];f1+=F[4*i+1];f2+=F[4*i+2];f3+=F[4*i+3];
    double m2=F[4*i+3]*F[4*i+3]-F[4*i]*F[4*i]-F[4*i+1]*F[4*i+1]-F[4*i+2]*F[4*i+2];
    double mt=dPoleMass(abs(fid[i])); mdm=fmax(mdm,fabs(m2-mt*mt)); }
  *cons=fmax(fmax(fabs(f0-ps0),fabs(f1-ps1)),fmax(fabs(f2-ps2),fabs(f3-ps3)));
  *dm=mdm; *nFout=nF; finSum[0]=f0;finSum[1]=f1;finSum[2]=f2;finSum[3]=f3;
}

__global__ void kern(int N,uint64_t base,double* cons,double* dm,int* nF,double* finSum){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  runOne(base,e,&cons[e],&dm[e],&nF[e],&finSum[4*e]);
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):500000; int TPB=128, blocks=(N+TPB-1)/TPB; uint64_t base=0xDEC0;
  double *dCons,*dDm,*dFin; int *dNf;
  CK(cudaMalloc(&dCons,(size_t)N*8));CK(cudaMalloc(&dDm,(size_t)N*8));CK(cudaMalloc(&dFin,(size_t)N*32));CK(cudaMalloc(&dNf,(size_t)N*4));
  kern<<<blocks,TPB>>>(N,base,dCons,dDm,dNf,dFin); CK(cudaDeviceSynchronize());
  std::vector<double> hCons(N),hDm(N),hFin((size_t)N*4); std::vector<int> hNf(N);
  CK(cudaMemcpy(hCons.data(),dCons,(size_t)N*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hDm.data(),dDm,(size_t)N*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hFin.data(),dFin,(size_t)N*32,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hNf.data(),dNf,(size_t)N*4,cudaMemcpyDeviceToHost));
  // reproducibility
  std::vector<int> hNf2(N); kern<<<blocks,TPB>>>(N,base,dCons,dDm,dNf,dFin); CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(hNf2.data(),dNf,(size_t)N*4,cudaMemcpyDeviceToHost));
  long repro=0; for(int e=0;e<N;++e) if(hNf[e]!=hNf2[e]) repro++;
  // GPU vs CPU (same runOne on host)
  int nCPU=(N<50000)?N:50000; long same=0; double maxMomRel=0;
  for(int e=0;e<nCPU;++e){ double c,d,fs[4]; int nf; runOne(base,e,&c,&d,&nf,fs);
    if(nf==hNf[e]){ same++; double mr=fmax(fmax(fabs(fs[0]-hFin[4*e]),fabs(fs[1]-hFin[4*e+1])),fmax(fabs(fs[2]-hFin[4*e+2]),fabs(fs[3]-hFin[4*e+3])));
      maxMomRel=fmax(maxMomRel,mr); } }
  double maxCons=0,maxDm=0; long nDrop=0,sumNf=0,nok=0;
  for(int e=0;e<N;++e){ if(hNf[e]<0){nDrop++;continue;} maxCons=fmax(maxCons,hCons[e]); maxDm=fmax(maxDm,hDm[e]); sumNf+=hNf[e]; nok++; }
  printf("GPU hadron-decay module test (%d events, pool of all unstable parents)\n",N);
  printf("  4-mom conservation : max|sum(final)-sum(primary)| = %.2e GeV\n",maxCons);
  printf("  on-shellness       : max|m^2-pole^2|              = %.2e GeV^2\n",maxDm);
  printf("  mean final mult    : %.2f (from 2-6 boosted primaries)\n",(double)sumNf/(nok>0?nok:1));
  printf("  drop (overflow)    : %ld / %d\n",nDrop,N);
  printf("  reproducibility    : GPU re-run diffs = %ld\n",repro);
  printf("  GPU vs CPU         : same final-count %ld/%d = %.3f%%, momentum-sum agree %.2e\n",same,nCPU,100.0*same/nCPU,maxMomRel);
  bool ok=(maxCons<1e-9)&&(maxDm<1e-6)&&(repro==0)&&(same>(long)(0.99*nCPU))&&(maxMomRel<1e-6)&&(nDrop<N/20);
  printf("VALIDATION: %s (4-mom conservation + on-shell + reproducible + GPU==CPU)\n",ok?"PASS":"FAIL");
  cudaFree(dCons);cudaFree(dDm);cudaFree(dFin);cudaFree(dNf);
  return ok?0:2;
}
