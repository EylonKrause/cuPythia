// cuPythia kernel 13 — O(N^2) all-pairs hadronic rescattering (heavy-ion).
//
// The Pythia core team flagged hadronic rescattering as the prime GPU target: in
// heavy-ion events thousands of hadrons are produced and the collision-finding
// cost grows as the SQUARE of the multiplicity (all hadron pairs are screened).
// That O(N^2) all-pairs screen is exactly the embarrassingly-parallel shape a GPU
// wants. Here we generate N hadrons at freeze-out and count interacting pairs
// (closest approach within the interaction radius d = sqrt(sigma/pi)) — on the GPU
// (one thread per hadron, inner loop over partners) and validate the count is
// IDENTICAL to a CPU reference.
//
// (This is the geometric collision-FINDING core. Processing collisions in time
// order is intrinsically sequential — the genuinely hard research piece — and is
// the next step; see README.)
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o rescatter rescatter.cu
// Run:   ./rescatter [nHadrons=20000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

// One thread per hadron i; count partners j>i within the interaction radius.
__global__ void pairKernel(const double* x,const double* y,const double* z,int N,
                           double d2max, unsigned long long* count){
  int i = blockIdx.x*blockDim.x + threadIdx.x; if(i>=N) return;
  double xi=x[i], yi=y[i], zi=z[i];
  unsigned long long c=0;
  for(int j=i+1;j<N;++j){
    double dx=xi-x[j], dy=yi-y[j], dz=zi-z[j];
    if(dx*dx+dy*dy+dz*dz < d2max) ++c;
  }
  atomicAdd(count,c);
}
static unsigned long long cpuCount(const std::vector<double>&x,const std::vector<double>&y,
                                   const std::vector<double>&z,int N,double d2max){
  unsigned long long c=0;
  for(int i=0;i<N;++i) for(int j=i+1;j<N;++j){
    double dx=x[i]-x[j],dy=y[i]-y[j],dz=z[i]-z[j];
    if(dx*dx+dy*dy+dz*dz < d2max) ++c;
  } return c;
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)
using Clock=std::chrono::steady_clock;

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):20000;
  double R=7.0;                       // fm, Pb-nucleus-scale fireball
  double sigma_fm2=4.0;               // 40 mb total cross section = 4 fm^2
  double d2max=sigma_fm2/M_PI;        // interaction radius^2 = sigma/pi
  // Deterministic freeze-out positions (same on host + device via SplitMix64).
  std::vector<double> x(N),y(N),z(N);
  uint64_t ctr=0x4EUL;
  for(int i=0;i<N;++i){ x[i]=(2.0*u01(splitmix64(ctr++))-1.0)*R;
                        y[i]=(2.0*u01(splitmix64(ctr++))-1.0)*R;
                        z[i]=(2.0*u01(splitmix64(ctr++))-1.0)*R; }

  double *dx,*dy,*dz; unsigned long long* dCount;
  CK(cudaMalloc(&dx,N*sizeof(double))); CK(cudaMalloc(&dy,N*sizeof(double))); CK(cudaMalloc(&dz,N*sizeof(double)));
  CK(cudaMalloc(&dCount,sizeof(unsigned long long))); CK(cudaMemset(dCount,0,sizeof(unsigned long long)));
  CK(cudaMemcpy(dx,x.data(),N*sizeof(double),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dy,y.data(),N*sizeof(double),cudaMemcpyHostToDevice));
  CK(cudaMemcpy(dz,z.data(),N*sizeof(double),cudaMemcpyHostToDevice));
  int threads=256, blocks=(N+threads-1)/threads;
  pairKernel<<<blocks,threads>>>(dx,dy,dz,N,d2max,dCount); CK(cudaDeviceSynchronize()); // warmup
  CK(cudaMemset(dCount,0,sizeof(unsigned long long)));
  auto g0=Clock::now();
  pairKernel<<<blocks,threads>>>(dx,dy,dz,N,d2max,dCount);
  CK(cudaDeviceSynchronize());
  auto g1=Clock::now();
  unsigned long long gpu=0; CK(cudaMemcpy(&gpu,dCount,sizeof(unsigned long long),cudaMemcpyDeviceToHost));

  auto c0=Clock::now(); unsigned long long cpu=cpuCount(x,y,z,N,d2max); auto c1=Clock::now();
  double gms=std::chrono::duration<double,std::milli>(g1-g0).count();
  double cms=std::chrono::duration<double,std::milli>(c1-c0).count();
  double pairs=0.5*(double)N*(N-1);

  printf("O(N^2) hadronic rescattering screen (heavy-ion), N=%d hadrons\n", N);
  printf("  pairs screened   = %.3e   interaction radius = %.3f fm\n", pairs, sqrt(d2max));
  printf("  interacting pairs: GPU = %llu   CPU = %llu\n", (unsigned long long)gpu,(unsigned long long)cpu);
  printf("  GPU %.1f ms (%.3e pairs/s)   CPU %.1f ms   speedup %.1fx\n",
         gms, pairs/(gms/1000.0), cms, cms/gms);
  bool ok = (gpu==cpu);
  printf("VALIDATION: %s (GPU all-pairs count == CPU reference, exact)\n", ok?"PASS":"FAIL");
  cudaFree(dx);cudaFree(dy);cudaFree(dz);cudaFree(dCount);
  return ok?0:2;
}
