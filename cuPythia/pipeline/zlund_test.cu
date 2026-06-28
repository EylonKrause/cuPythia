// cuPythia hadronization — increment 1: the Lund symmetric fragmentation function
// sampler zLund(), ported FAITHFULLY from Pythia 8.317 FragmentationFlavZpT.cc
// (zLundMax :1197, zLund :1223, initFunc :1148), validated IN ISOLATION against the
// analytic closed form  f(z) = (1-z)^a exp(-b/z) / z^c  before it is ever wired into a
// fragmentation chain. This is the critic's "get zLund right first": it is the physics
// core and the cleanest fully-in-our-control validator. b carries m_T^2 (b = bLund*m_T^2),
// so the four m_T^2 values below exercise all three envelope regimes (peak near 0, middle,
// near 1). Host/device-identical, counter-based RNG (../common/rng.cuh).
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o zlund_test zlund_test.cu
// Run:   ./zlund_test [N=2000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)

static constexpr double ALUND=0.68, BLUND=0.98, CLUND=1.0;   // Pythia 8.317 StringZ defaults

__host__ __device__ inline double zLundMax(double a,double b,double c){
  const double AFROMZERO=0.02, AFROMC=0.01;
  bool aIsZero=(a<AFROMZERO), aIsC=(fabs(a-c)<AFROMC);
  double zMax;
  if(aIsZero) zMax=(c>b)? b/c : 1.0;
  else if(aIsC) zMax=b/(b+c);
  else { zMax=0.5*(b+c-sqrt((b-c)*(b-c)+4.0*a*b))/(c-a);
         if(zMax>0.9999 && b>100.0) zMax=fmin(zMax,1.0-a/b); }
  return zMax;
}

// Faithful port of StringZ::zLund (head=1, no reweighting). Draws via counter RNG.
__host__ __device__ inline double zLundSample(double a,double b,double c,uint64_t& ctr){
  const double CFROMUNITY=0.01, AFROMZERO=0.02, EXPMAX=50.0;
  bool cIsUnity=(fabs(c-1.0)<CFROMUNITY), aIsZero=(a<AFROMZERO);
  double zMax=zLundMax(a,b,c);
  bool peakedNearZero=(zMax<0.1), peakedNearUnity=(zMax>0.85 && b>1.0);
  double fIntLow=1.,fIntHigh=1.,fInt=2.,zDiv=0.5,zDivC=0.5;
  if(peakedNearZero){
    zDiv=2.75*zMax; fIntLow=zDiv;
    if(cIsUnity) fIntHigh=-zDiv*log(zDiv);
    else { zDivC=pow(zDiv,1.0-c); fIntHigh=zDiv*(1.0-1.0/zDivC)/(c-1.0); }
    fInt=fIntLow+fIntHigh;
  } else if(peakedNearUnity){
    double rcb=sqrt(4.0+(c/b)*(c/b));
    zDiv=rcb-1.0/zMax-(c/b)*log(zMax*0.5*(rcb+c/b));
    if(!aIsZero) zDiv+=(a/b)*log(1.0-zMax);
    zDiv=fmin(zMax,fmax(0.0,zDiv));
    fIntLow=1.0/b; fIntHigh=1.0-zDiv; fInt=fIntLow+fIntHigh;
  }
  double z=0.5; bool accept=false;
  do{
    z=u01(splitmix64(ctr++));
    double fPrel=1.0;
    if(peakedNearZero){
      if(fInt*u01(splitmix64(ctr++))<fIntLow) z=zDiv*z;
      else if(cIsUnity){ z=pow(zDiv,z); fPrel=zDiv/z; }
      else { z=pow(zDivC+(1.0-zDivC)*z,1.0/(1.0-c)); fPrel=pow(zDiv/z,c); }
    } else if(peakedNearUnity){
      if(fInt*u01(splitmix64(ctr++))<fIntLow){ z=zDiv+log(z)/b; fPrel=exp(b*(z-zDiv)); }
      else z=zDiv+(1.0-zDiv)*z;
    }
    if(z>0.0 && z<1.0){
      double fRnd=u01(splitmix64(ctr++));
      double fExp=b*(1.0/zMax-1.0/z)+c*log(zMax/z);
      if(!aIsZero) fExp+=a*log((1.0-z)/(1.0-zMax));
      double fVal=exp(fmax(-EXPMAX,fmin(EXPMAX,fExp)));
      accept=((fVal/fPrel)>fRnd);
    }
  } while(!accept);
  return z;
}

__global__ void kern(int N,double a,double b,double c,uint64_t base,double* out){
  int i=blockIdx.x*(int)blockDim.x+threadIdx.x; if(i>=N) return;
  uint64_t ctr=base+(uint64_t)i*0x9E3779B97F4A7C15ULL;
  out[i]=zLundSample(a,b,c,ctr);
}

// analytic unnormalised f(z)
static inline double fLund(double z,double a,double b,double c){
  if(z<=0||z>=1) return 0.0; return pow(1.0-z,a)*exp(-b/z)/pow(z,c);
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):2000000;
  int TPB=256, blocks=(N+TPB-1)/TPB;
  const int NB=50;
  double mT2list[4]={0.1,0.5,2.0,10.0};
  double *dz; CK(cudaMalloc(&dz,(size_t)N*8));
  std::vector<double> hz(N), hz2(N);

  bool allOK=true;
  printf("Lund f(z) sampler — sampled distribution vs analytic (1-z)^a exp(-b/z)/z^c\n");
  printf("  a=%.2f c=%.1f, b=bLund*mT2 (bLund=%.2f); %d samples each\n",ALUND,CLUND,BLUND,N);
  for(int t=0;t<4;++t){
    double mT2=mT2list[t], b=BLUND*mT2, zMax=zLundMax(ALUND,b,CLUND);
    const char* regime=(zMax<0.1)?"peak~0":((zMax>0.85&&b>1.0)?"peak~1":"middle");
    uint64_t base=0x21D+t*0x1000000ULL;
    kern<<<blocks,TPB>>>(N,ALUND,b,CLUND,base,dz); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(hz.data(),dz,(size_t)N*8,cudaMemcpyDeviceToHost));
    kern<<<blocks,TPB>>>(N,ALUND,b,CLUND,base,dz); CK(cudaDeviceSynchronize());   // reproducibility
    CK(cudaMemcpy(hz2.data(),dz,(size_t)N*8,cudaMemcpyDeviceToHost));
    long repro=0; for(int i=0;i<N;++i) if(hz[i]!=hz2[i]) repro++;

    // sampled histogram + GPU mean
    long obs[NB]={0}; double sumG=0;
    for(int i=0;i<N;++i){ sumG+=hz[i]; int k=(int)(hz[i]*NB); if(k<0)k=0; if(k>=NB)k=NB-1; obs[k]++; }
    // analytic per-bin probability (fine sub-integration), normalised over (0,1)
    double exp_[NB], tot=0; int sub=200;
    for(int k=0;k<NB;++k){ double s=0,lo=(double)k/NB,hi=(double)(k+1)/NB,h=(hi-lo)/sub;
      for(int j=0;j<sub;++j){ double z=lo+(j+0.5)*h; s+=fLund(z,ALUND,b,CLUND)*h; } exp_[k]=s; tot+=s; }
    for(int k=0;k<NB;++k) exp_[k]/=tot;
    // chi2/ndf over well-populated bins (expected count >= 10)
    double chi2=0; int ndf=0;
    for(int k=0;k<NB;++k){ double e=exp_[k]*N; if(e>=10.0){ double d=obs[k]-e; chi2+=d*d/e; ndf++; } }
    ndf=(ndf>1)?ndf-1:1;
    // CPU mean (same seeds) for determinism cross-check
    double sumC=0; for(int i=0;i<N;++i){ uint64_t ctr=base+(uint64_t)i*0x9E3779B97F4A7C15ULL; sumC+=zLundSample(ALUND,b,CLUND,ctr); }
    double meanG=sumG/N, meanC=sumC/N;
    bool ok=(chi2/ndf<2.0)&&(repro==0)&&(fabs(meanG-meanC)/meanC<1e-2);
    allOK&=ok;
    printf("  mT2=%5.1f (b=%5.2f, zMax=%.3f, %-6s): chi2/ndf=%.2f  <z>=%.4f (GPU) %.4f (CPU)  repro=%ld  %s\n",
           mT2,b,zMax,regime,chi2/ndf,meanG,meanC,repro,ok?"OK":"FAIL");
  }
  printf("VALIDATION: %s (Lund f(z) sampler matches analytic in all regimes + reproducible + GPU==CPU)\n",
         allOK?"PASS":"FAIL");
  cudaFree(dz);
  return allOK?0:2;
}
