// cuPythia kernel 10 — on-GPU event compaction with CUB (scalable I/O).
//
// A concrete example of a CUDA toolkit library helping cuPythia. Kernel 07
// filtered accepted (unweighted) events on the HOST. Here we generate candidate
// events on the GPU, flag the accepted ones, and use cub::DeviceSelect::Flagged
// to compact them into a dense array ENTIRELY ON THE GPU — no host round-trip to
// filter. This is the scalable unweighted-event I/O pattern modern GPU generators
// (madgraph4gpu) rely on. CUB ships with the CUDA toolkit (header-only), so this
// builds on Windows and Linux with no extra dependency.
//
// Build: nvcc -O3 -arch=sm_120 -o cub_compaction cub_compaction.cu
// Run:   ./cub_compaction [nCandidates=16777216]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include "../common/rng.cuh"

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

// One thread per candidate: sample cosθ, accept with prob w/wMax, emit flag.
__global__ void genKernel(uint64_t seed,int n,double s,double alpS,double cMax,double wMax,
                          double* cosOut, unsigned char* flags, unsigned long long* accCount){
  int i = blockIdx.x*blockDim.x + threadIdx.x; if(i>=n) return;
  uint64_t ctr = seed + (uint64_t)i*0x100000001B3ULL;
  double c=(2.0*u01(splitmix64(ctr++))-1.0)*cMax;
  double w=weightAt(c,s,alpS);
  double u=u01(splitmix64(ctr++));
  unsigned char f=(u*wMax < w)?1:0;
  cosOut[i]=c; flags[i]=f;
  if(f) atomicAdd(accCount,1ULL);
}

static double wMaxScan(int N,double s,double alpS,double cMax){
  double m=0; for(int i=0;i<=N;++i){ double c=-cMax+(2.0*cMax)*i/N; m=fmax(m,weightAt(c,s,alpS)); } return m;
}
static double simpson(int N,double s,double alpS,double cMax){
  double a=-cMax,h=(2.0*cMax)/N,sum=0;
  for(int i=0;i<=N;++i){ double c=a+i*h,w=(i==0||i==N)?1.0:(i%2?4.0:2.0); sum+=w*weightAt(c,s,alpS); }
  return (sum*h/3.0)*(s/2.0);
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  int n = (argc>1)? atoi(argv[1]) : (1<<24);   // 16.7M candidates
  double s=100.0*100.0, alphaS=0.118, cMax=0.9, conv_pb=0.3893793721e9;
  uint64_t seed=0xC0B0ULL;
  double wMax=wMaxScan(2000000,s,alphaS,cMax);

  double* dCos; unsigned char* dFlags; double* dOut; unsigned long long* dAcc; int* dNum;
  CK(cudaMalloc(&dCos,(size_t)n*sizeof(double)));
  CK(cudaMalloc(&dFlags,(size_t)n*sizeof(unsigned char)));
  CK(cudaMalloc(&dOut,(size_t)n*sizeof(double)));
  CK(cudaMalloc(&dAcc,sizeof(unsigned long long))); CK(cudaMemset(dAcc,0,sizeof(unsigned long long)));
  CK(cudaMalloc(&dNum,sizeof(int)));
  int threads=256, blocks=(n+threads-1)/threads;
  genKernel<<<blocks,threads>>>(seed,n,s,alphaS,cMax,wMax,dCos,dFlags,dAcc);
  CK(cudaDeviceSynchronize());

  // CUB stream compaction: keep dCos[i] where dFlags[i] != 0  ->  dOut, count -> dNum.
  void* dTemp=nullptr; size_t tempBytes=0;
  cub::DeviceSelect::Flagged(dTemp,tempBytes,dCos,dFlags,dOut,dNum,n);
  CK(cudaMalloc(&dTemp,tempBytes));
  cub::DeviceSelect::Flagged(dTemp,tempBytes,dCos,dFlags,dOut,dNum,n);
  CK(cudaDeviceSynchronize());

  int hNum=0; unsigned long long hAcc=0;
  CK(cudaMemcpy(&hNum,dNum,sizeof(int),cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(&hAcc,dAcc,sizeof(unsigned long long),cudaMemcpyDeviceToHost));

  double eta=(double)hNum/(double)n;
  double V=(2.0*cMax)*(s/2.0);
  double sigma_pb = V*eta*wMax*conv_pb;             // unweighted-count cross section
  double ref_pb = simpson(2000000,s,alphaS,cMax)*conv_pb;

  printf("On-GPU event compaction with CUB (cub::DeviceSelect::Flagged)\n");
  printf("  candidates generated      = %d\n", n);
  printf("  accepted (CUB-compacted)  = %d   (independent atomic count = %llu)\n",
         hNum, (unsigned long long)hAcc);
  printf("  unweighting efficiency    = %.2f%%\n", 100.0*eta);
  printf("  sigma from compacted set  = %.6e pb   Simpson ref = %.6e pb   relerr = %.2e\n",
         sigma_pb, ref_pb, fabs(sigma_pb-ref_pb)/ref_pb);
  bool ok = ((unsigned long long)hNum==hAcc) && fabs(sigma_pb-ref_pb)/ref_pb < 3e-3;
  printf("VALIDATION: %s (CUB count == independent count AND sigma matches quadrature)\n",
         ok?"PASS":"FAIL");
  cudaFree(dCos);cudaFree(dFlags);cudaFree(dOut);cudaFree(dAcc);cudaFree(dNum);cudaFree(dTemp);
  return ok?0:2;
}
