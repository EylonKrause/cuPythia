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

#include "shower_inc.cuh"

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
#ifdef ME_FIRST
  bool detOK = (structSame > (long)(0.999*nCPU));  // ME adds a transcendental veto -> rare ULP flips OK (decided a priori)
#else
  bool detOK = (structSame == nCPU);               // pure LL: exact bit-identical control flow
#endif
  bool ok = (maxMomViol<1e-5) && (maxM2<1e-3) && (reproDiff==0) &&
            (meanN>2.5&&meanN<25.0) && detOK && (maxMomRel<1e-6) &&
            (mean1mT>0.01 && mean1mT<0.30);
  printf("VALIDATION: %s (momentum+on-shell+reproducible+CPU-agreement+thrust-sane)\n", ok?"PASS":"FAIL");
  cudaFree(dN);cudaFree(dTot);cudaFree(dM2);cudaFree(dThr);
  return ok?0:2;
}
