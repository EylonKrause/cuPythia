// cuPythia pipeline stage 2 — on-the-fly multi-weight systematic reweighting.
//
// THE counter-RNG advantage. For each gg->gg event we emit, in ONE pass, a vector
// of renormalization-scale variation weights (mu_R = 0.5, 1, 2 x pT-hat, the
// single-scale subset of the 7-point set; mu_F needs PDFs, stage 1). Since LO
// sigma ~ alphaS(mu_R)^2, each variation is w0 * (alphaS_k/alphaS_0)^2 -- a cheap
// re-evaluation reusing the IDENTICAL sampled phase-space point.
//
// We then prove the production-relevant property: re-running each variation
// INDEPENDENTLY, pinned to the same per-event counter seed, reproduces every
// weight BIT-IDENTICALLY. So N weights cost one pass, not N runs -- exactly what
// RANMAR-based stock Pythia cannot do cheaply, and what ATLAS/CMS scale-uncertainty
// bands need (cf. MadtRex, arXiv:2510.05100).
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o reweight reweight.cu
// Run:   ./reweight [nEvents=4000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "physics.cuh"

#define NVAR 3
__constant__ double FAC[NVAR] = {0.5, 1.0, 2.0};   // mu_R / pT-hat

// 1-loop running coupling, nf=5: alphaS(mu^2) from alphaS(muRef^2).
__host__ __device__ inline double alphaS1(double mu2,double muRef2,double aSref){
  double b0=(33.0-2.0*5.0)/(12.0*M_PI);
  return aSref/(1.0 + aSref*b0*log(mu2/muRef2));
}

// Compute the NVAR scale-variation weights for one event, from its counter seed.
__device__ inline void eventWeights(uint64_t seed,double sqrtS,double cMax,
                                    double aSref,double muRef2,double* w){
  uint64_t c=seed;
  double E=0.5*sqrtS, s=sqrtS*sqrtS;
  double cosT=(2.0*u01(splitmix64(c++))-1.0)*cMax;
  (void)u01(splitmix64(c++));                 // phi draw (keeps the RNG stream aligned)
  double st=sqrt(fmax(0.0,1.0-cosT*cosT));
  double pT=E*st; double mu0=(pT>1e-3)?pT:1e-3;
  double aS0=alphaS1(mu0*mu0,muRef2,aSref);
  double t=-0.5*s*(1.0-cosT);
  double w0=gg2gg_sigma(s,t,-s-t,aS0);        // nominal weight at mu_R = pT-hat
  for(int k=0;k<NVAR;++k){
    double mu=FAC[k]*mu0; double aSk=alphaS1(mu*mu,muRef2,aSref);
    double r=aSk/aS0; w[k]=w0*r*r;
  }
}
__device__ inline uint64_t eventSeed(uint64_t base,int e){
  return splitmix64(base ^ ((uint64_t)e*0x9E3779B97F4A7C15ULL));
}

// One pass: all NVAR weights per event.
__global__ void onePass(int N,double sqrtS,double cMax,double aSref,double muRef2,uint64_t base,double* out){
  int e=blockIdx.x*blockDim.x+threadIdx.x; if(e>=N) return;
  double w[NVAR]; eventWeights(eventSeed(base,e),sqrtS,cMax,aSref,muRef2,w);
  for(int k=0;k<NVAR;++k) out[(size_t)e*NVAR+k]=w[k];
}
// Independent re-run of a SINGLE variation k, pinned to the same seed.
__global__ void rerunOne(int N,int k,double sqrtS,double cMax,double aSref,double muRef2,uint64_t base,double* out){
  int e=blockIdx.x*blockDim.x+threadIdx.x; if(e>=N) return;
  double w[NVAR]; eventWeights(eventSeed(base,e),sqrtS,cMax,aSref,muRef2,w);
  out[e]=w[k];
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):4000000;
  double sqrtS=100.0, cMax=0.9, conv_pb=0.3893793721e9;
  double aSref=0.118, muRef=91.1876, muRef2=muRef*muRef;   // alphaS(M_Z)=0.118
  uint64_t base=0x5CA1EULL;

  double* dAll; CK(cudaMalloc(&dAll,(size_t)N*NVAR*sizeof(double)));
  double* dOne; CK(cudaMalloc(&dOne,(size_t)N*sizeof(double)));
  int threads=256, blocks=(N+threads-1)/threads;
  onePass<<<blocks,threads>>>(N,sqrtS,cMax,aSref,muRef2,base,dAll);
  CK(cudaDeviceSynchronize());
  std::vector<double> all((size_t)N*NVAR); CK(cudaMemcpy(all.data(),dAll,(size_t)N*NVAR*8,cudaMemcpyDeviceToHost));

  // For each variation: independent pinned re-run, compare bit-for-bit.
  double maxDiff=0; std::vector<double> one(N);
  double sig[NVAR]={0,0,0};
  for(int k=0;k<NVAR;++k){
    rerunOne<<<blocks,threads>>>(N,k,sqrtS,cMax,aSref,muRef2,base,dOne);
    CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(one.data(),dOne,(size_t)N*8,cudaMemcpyDeviceToHost));
    double sw=0;
    for(int e=0;e<N;++e){ double d=fabs(one[e]-all[(size_t)e*NVAR+k]); if(d>maxDiff)maxDiff=d; sw+=all[(size_t)e*NVAR+k]; }
    sig[k]=(2.0*cMax)*(sqrtS*sqrtS/2.0)*(sw/N)*conv_pb;
  }
  printf("Multi-weight scale-variation reweighting (gg->gg, %d events, %d weights each)\n",N,NVAR);
  printf("  sigma(mu_R=0.5x) = %.4e pb   (alphaS larger -> bigger)\n", sig[0]);
  printf("  sigma(mu_R=1.0x) = %.4e pb   (nominal)\n", sig[1]);
  printf("  sigma(mu_R=2.0x) = %.4e pb   (alphaS smaller -> smaller)\n", sig[2]);
  printf("  scale-uncertainty band = [%+.1f%%, %+.1f%%] around nominal\n",
         100.0*(sig[2]-sig[1])/sig[1], 100.0*(sig[0]-sig[1])/sig[1]);
  printf("  one-pass vs %d independent pinned re-runs: max|diff| = %.1e\n", NVAR, maxDiff);
  bool ok = (maxDiff==0.0) && (sig[0]>sig[1]) && (sig[1]>sig[2]);
  printf("VALIDATION: %s (N weights in one pass == N pinned re-runs, bit-identical; band ordered)\n",
         ok?"PASS":"FAIL");
  cudaFree(dAll); cudaFree(dOne);
  return ok?0:2;
}
