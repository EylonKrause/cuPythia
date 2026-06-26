// cuPythia kernel 14 — batched parton shower (Sudakov veto), the hard one.
//
// The parton shower is the SEQUENTIAL heart of event generation: emissions happen
// in order, each depending on the last, so it does not vectorise within one event.
// The GPU pattern is to batch ACROSS events: one shower per thread, thousands of
// independent showers in flight. Each shower evolves a scale t = pT^2 downward
// from t_max via the Sudakov form factor; for a 1/t kernel with constant
// integrated splitting C, the next scale is the exact inversion t -> t * R^(1/C).
//
// The clean validation is the SUDAKOV no-emission probability between t_max and
// t_min:   Delta = exp(-C * ln(t_max/t_min)),   and the emission multiplicity is
// Poisson with mean  C * ln(t_max/t_min).  We reproduce both on GPU.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o shower shower.cu
// Run:   ./shower [nShowers=8000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

// One thread = one shower. Evolve t downward by the Sudakov inversion; count
// emissions above t_min.  (Sequential within a shower; batched across the grid.)
__global__ void showerKernel(uint64_t seed,int nShowers,double C,double tMax,double tMin,
                             unsigned long long* zeroCount, unsigned long long* totEmit){
  int s = blockIdx.x*blockDim.x + threadIdx.x; if(s>=nShowers) return;
  uint64_t ctr = seed + (uint64_t)s*0x100000001B3ULL;
  double t=tMax; unsigned n=0;
  for(;;){
    double R = u01(splitmix64(ctr++));
    t *= pow(R, 1.0/C);          // next emission scale (Sudakov form-factor inversion)
    if(t < tMin) break;          // dropped below the cutoff -> shower terminates
    ++n;
  }
  if(n==0) atomicAdd(zeroCount,1ULL);
  atomicAdd(totEmit,(unsigned long long)n);
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  int nShowers=(argc>1)?atoi(argv[1]):8000000;
  double alphaS=0.118, zcut=0.1;
  // Effective integrated quark splitting (4/3)*int_{zcut}^{1-zcut} (1+z^2)/(1-z) dz
  // ~ (4/3)*(2 ln(1/zcut) - 3/2);  C = (alphaS/2pi) * that.
  double splitInt = (4.0/3.0)*(2.0*log(1.0/zcut) - 1.5);
  double C = (alphaS/(2.0*M_PI))*splitInt;
  double tMax=1.0e4, tMin=1.0, L=log(tMax/tMin);   // pT^2 in GeV^2

  unsigned long long *dZero,*dTot;
  CK(cudaMalloc(&dZero,8)); CK(cudaMalloc(&dTot,8)); CK(cudaMemset(dZero,0,8)); CK(cudaMemset(dTot,0,8));
  int threads=256, blocks=(nShowers+threads-1)/threads;
  showerKernel<<<blocks,threads>>>(0x540ULL,nShowers,C,tMax,tMin,dZero,dTot);
  CK(cudaDeviceSynchronize());
  unsigned long long zero=0,tot=0;
  CK(cudaMemcpy(&zero,dZero,8,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(&tot,dTot,8,cudaMemcpyDeviceToHost));

  double fracZero=(double)zero/nShowers, meanN=(double)tot/nShowers;
  double sudakov=exp(-C*L), meanAna=C*L;
  printf("Batched parton shower (Sudakov veto), %d showers\n", nShowers);
  printf("  C = %.4f   ln(tMax/tMin) = %.3f   (tMax=%.0f, tMin=%.0f GeV^2)\n", C, L, tMax, tMin);
  printf("  no-emission fraction: MC = %.5f   Sudakov exp(-C L) = %.5f   (rel %.2e)\n",
         fracZero, sudakov, fabs(fracZero-sudakov)/sudakov);
  printf("  mean emissions:       MC = %.5f   analytic C L      = %.5f   (rel %.2e)\n",
         meanN, meanAna, fabs(meanN-meanAna)/meanAna);
  bool ok = fabs(fracZero-sudakov)/sudakov < 5e-3 && fabs(meanN-meanAna)/meanAna < 5e-3;
  printf("VALIDATION: %s (no-emission prob == Sudakov AND mean mult == C ln)\n", ok?"PASS":"FAIL");
  cudaFree(dZero);cudaFree(dTot);
  return ok?0:2;
}
