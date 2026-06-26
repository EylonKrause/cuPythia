// cuPythia kernel 11 — VEGAS adaptive importance sampling (raises unweighting eff).
//
// Unweighting efficiency eta = <w>/w_max is the cost driver: kernel 07 measured
// eta ~ 10% for gg->gg with UNIFORM sampling, i.e. 90% of generated events are
// thrown away. VEGAS (Lepage 1978) adapts a piecewise sampling grid to the
// integrand so the MC weight w = f/p flattens, pushing eta up. This is the
// classical precursor to neural importance sampling (MadNIS, arXiv:2212.06172).
//
// Per iteration: a GPU kernel samples from the current grid, accumulates the
// integral, the max weight, and per-bin importance (shared-memory binning); the
// host rebins the grid for equal importance. eta is reported before vs after.
//
// Build: nvcc -O3 -arch=sm_120 -o vegas vegas.cu
// Run:   ./vegas [samplesPerThread=2000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

#define K 128   // VEGAS grid bins

__host__ __device__ inline double pow2(double x){ return x*x; }
__host__ __device__ inline double gg2gg_sigma(double sH,double tH,double uH,double alpS){
  double s2=sH*sH,t2=tH*tH,u2=uH*uH;
  double a=(9./4.)*(t2/s2+2.*tH/sH+3.+2.*sH/tH+s2/t2);
  double b=(9./4.)*(u2/s2+2.*uH/sH+3.+2.*sH/uH+s2/u2);
  double c=(9./4.)*(t2/u2+2.*tH/uH+3.+2.*uH/tH+u2/t2);
  return (M_PI/s2)*pow2(alpS)*0.5*(a+b+c);
}
__host__ __device__ inline double weightAt(double cosT,double s,double alpS){
  double t=-0.5*s*(1.0-cosT); return gg2gg_sigma(s,t,-s-t,alpS);
}

// Sample from the current grid; accumulate integral, max weight, per-bin importance.
__global__ void vegasKernel(uint64_t seed, uint64_t Tper, const double* grid,
                            double s, double alpS,
                            double* gSumW, double* gSumW2,
                            unsigned long long* gMaxW, double* gDacc){
  __shared__ double sdacc[K];
  for(int t=threadIdx.x;t<K;t+=blockDim.x) sdacc[t]=0.0;
  __syncthreads();
  uint64_t tid=blockIdx.x*(uint64_t)blockDim.x+threadIdx.x;
  uint64_t ctr=seed + tid*0x100000001B3ULL;
  double lsum=0.0,lsum2=0.0; unsigned long long lmax=0ULL;
  for(uint64_t i=0;i<Tper;++i){
    double xi=u01(splitmix64(ctr++));
    int b=(int)(xi*K); if(b>=K) b=K-1;
    double lo=grid[b], hi=grid[b+1], d=hi-lo;
    double c=lo + u01(splitmix64(ctr++))*d;
    double f=weightAt(c,s,alpS);
    double w=f*(double)K*d;                 // w = f / pdf,  pdf = (1/K)/d
    lsum+=w; lsum2+=w*w;
    unsigned long long wb=(unsigned long long)__double_as_longlong(w);
    if(wb>lmax) lmax=wb;
    atomicAdd(&sdacc[b], f*d);               // bin's share of the integral
  }
  __syncthreads();
  atomicAdd(gSumW,lsum); atomicAdd(gSumW2,lsum2);
  atomicMax(gMaxW,lmax);
  for(int t=threadIdx.x;t<K;t+=blockDim.x) atomicAdd(&gDacc[t],sdacc[t]);
}

// VEGAS rebin: redistribute edges so each bin carries equal (smoothed) importance.
static void rebin(double* grid, const double* dacc){
  double d[K];
  for(int i=0;i<K;i++){
    double l=(i>0)?dacc[i-1]:dacc[i], r=(i<K-1)?dacc[i+1]:dacc[i];
    d[i]=(l+dacc[i]+r)/3.0 + 1e-300;        // light smoothing, avoid zeros
  }
  double tot=0; for(int i=0;i<K;i++) tot+=d[i];
  double step=tot/K, ng[K+1]; ng[0]=grid[0]; ng[K]=grid[K];
  double cur=0; int oj=0;
  for(int k=1;k<K;k++){
    double tgt=k*step;
    while(oj<K-1 && cur+d[oj]<tgt){ cur+=d[oj]; oj++; }
    double frac=(tgt-cur)/d[oj]; if(frac<0)frac=0; if(frac>1)frac=1;
    ng[k]=grid[oj]+frac*(grid[oj+1]-grid[oj]);
  }
  for(int i=0;i<=K;i++) grid[i]=ng[i];
}
static double simpson(int N,double s,double alpS,double cMax){
  double a=-cMax,h=(2.0*cMax)/N,sum=0;
  for(int i=0;i<=N;++i){ double c=a+i*h,w=(i==0||i==N)?1.0:(i%2?4.0:2.0); sum+=w*weightAt(c,s,alpS); }
  return (sum*h/3.0)*(s/2.0);
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  uint64_t Tper=(argc>1)?strtoull(argv[1],nullptr,10):2000ULL;
  const int blocks=256, threads=256, NITER=12;
  double s=100.0*100.0, alphaS=0.118, cMax=0.9, conv_pb=0.3893793721e9;
  double M=(double)blocks*threads*Tper;

  double grid[K+1]; for(int i=0;i<=K;i++) grid[i]=-cMax+(2.0*cMax)*i/K;
  double *dGrid,*dSumW,*dSumW2,*dDacc; unsigned long long* dMaxW;
  CK(cudaMalloc(&dGrid,(K+1)*sizeof(double)));
  CK(cudaMalloc(&dSumW,sizeof(double))); CK(cudaMalloc(&dSumW2,sizeof(double)));
  CK(cudaMalloc(&dDacc,K*sizeof(double))); CK(cudaMalloc(&dMaxW,sizeof(unsigned long long)));

  double etaFirst=0, etaLast=0, Ilast=0;
  for(int it=0; it<NITER; ++it){
    CK(cudaMemcpy(dGrid,grid,(K+1)*sizeof(double),cudaMemcpyHostToDevice));
    CK(cudaMemset(dSumW,0,sizeof(double))); CK(cudaMemset(dSumW2,0,sizeof(double)));
    CK(cudaMemset(dDacc,0,K*sizeof(double))); CK(cudaMemset(dMaxW,0,sizeof(unsigned long long)));
    vegasKernel<<<blocks,threads>>>(0xBEEF0000ULL + (uint64_t)it, Tper, dGrid, s, alphaS, dSumW, dSumW2, dMaxW, dDacc);
    CK(cudaDeviceSynchronize());
    double sumW=0,dacc[K]; unsigned long long maxWb=0;
    CK(cudaMemcpy(&sumW,dSumW,sizeof(double),cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&maxWb,dMaxW,sizeof(unsigned long long),cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(dacc,dDacc,K*sizeof(double),cudaMemcpyDeviceToHost));
    double meanW=sumW/M, maxW; memcpy(&maxW,&maxWb,sizeof(double));
    double eta=meanW/maxW;
    Ilast=meanW*(s/2.0)*conv_pb;   // <w> = integral d(cosθ); * (s/2) Jacobian -> sigma

    if(it==0) etaFirst=eta;
    etaLast=eta;
    rebin(grid,dacc);
  }
  double ref=simpson(2000000,s,alphaS,cMax)*conv_pb;

  printf("VEGAS adaptive importance sampling (gg->gg, %d bins, %d iters)\n",K,NITER);
  printf("  unweighting efficiency: uniform = %.2f%%  ->  VEGAS = %.2f%%   (%.1fx better)\n",
         100.0*etaFirst, 100.0*etaLast, etaLast/etaFirst);
  printf("  integral (VEGAS) = %.6e pb   Simpson ref = %.6e pb   relerr = %.2e\n",
         Ilast, ref, fabs(Ilast-ref)/ref);
  bool ok = (etaLast > 1.5*etaFirst) && (fabs(Ilast-ref)/ref < 5e-3);
  printf("VALIDATION: %s (VEGAS raises efficiency AND integral matches quadrature)\n", ok?"PASS":"FAIL");
  cudaFree(dGrid);cudaFree(dSumW);cudaFree(dSumW2);cudaFree(dDacc);cudaFree(dMaxW);
  return ok?0:2;
}
