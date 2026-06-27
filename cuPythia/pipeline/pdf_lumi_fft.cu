// cuPythia pipeline — FFT-accelerated PARTON LUMINOSITY (answer to: "can we FFT the convolution?").
//
// The hadronic cross section sigma = int dx1 dx2 dcos f(x1) f(x2) dsig(x1 x2 S) is a
// product-form MC integral, NOT a shift-convolution -> FFT does not apply to it, and
// an event generator must *sample* (x1,x2) per event anyway (a transform yields the
// integral, not events). BUT the parton luminosity
//     L(tau) = int_tau^1 (dx/x) f1(x) f2(tau/x)
// IS a multiplicative (Mellin) convolution. With x=e^-xi, tau=e^-eta it becomes an
// ADDITIVE convolution  L(eta) = int_0^eta F1(xi) F2(eta-xi) dxi,  F_i(xi)=f_i(e^-xi),
// which cuFFT evaluates in O(N log N) instead of the direct O(N^2). This file proves
// the FFT path reproduces the direct convolution bit-for-bit (to FP/FFT roundoff) and
// is dramatically faster at large grids. It is an INCLUSIVE tool (total luminosity /
// cross section), separate from the exclusive event generator.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o pdf_lumi_fft pdf_lumi_fft.cu -lcufft
// Run:   ./pdf_lumi_fft [log2N=15]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <cufft.h>

#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)
#define FK(call) do{ cufftResult r_=(call); if(r_!=CUFFT_SUCCESS){ \
  printf("cuFFT error %d at %s:%d\n",(int)r_,__FILE__,__LINE__); return 1; } }while(0)

// Toy single-parton density f(x)=4 x^0.5 (1-x)^3 -> F(xi)=f(e^-xi), zero-padded to length M.
__global__ void fillF(double* F,int N,int M,double h){
  int j=blockIdx.x*(int)blockDim.x+threadIdx.x; if(j>=M) return;
  if(j<N){ double x=exp(-(double)j*h); F[j]=4.0*pow(x,0.5)*pow(1.0-x,3.0); }
  else F[j]=0.0;
}
// Direct discrete linear convolution C[k] = sum_j F[j] F[k-j]  (O(N) per k -> O(N^2) total).
__global__ void directConv(const double* F,double* C,int N,int K){
  int k=blockIdx.x*(int)blockDim.x+threadIdx.x; if(k>=K) return;
  int jlo=(k-(N-1)>0)?k-(N-1):0, jhi=(k<N-1)?k:N-1; double s=0.0;
  for(int j=jlo;j<=jhi;++j) s+=F[j]*F[k-j];
  C[k]=s;
}
// Complex square A[i]^2 (F1==F2), in place into Cc.
__global__ void cSquare(const cufftDoubleComplex* A,cufftDoubleComplex* Cc,int n){
  int i=blockIdx.x*(int)blockDim.x+threadIdx.x; if(i>=n) return;
  double re=A[i].x, im=A[i].y; Cc[i].x=re*re-im*im; Cc[i].y=2.0*re*im;
}
__global__ void scaleR(double* x,int M,double s){ int i=blockIdx.x*(int)blockDim.x+threadIdx.x; if(i<M) x[i]*=s; }

int main(int argc,char**argv){
  int log2N=(argc>1)?atoi(argv[1]):15;
  int N=1<<log2N;                 // grid points on xi in [0,L]
  int K=2*N-1;                    // length of the linear convolution
  int M=1; while(M<K) M<<=1;      // FFT length (power of two, >= K) for zero-padded linear conv
  double L=12.0, h=L/(N-1);       // xi range; eta=xi grid carries tau=e^-eta down to e^-12 ~ 6e-6
  int Mc=M/2+1;                   // r2c output length
  int TPB=256;

  double *dF,*dCd,*dOut; CK(cudaMalloc(&dF,(size_t)M*8)); CK(cudaMalloc(&dCd,(size_t)K*8)); CK(cudaMalloc(&dOut,(size_t)M*8));
  cufftDoubleComplex *dA,*dC; CK(cudaMalloc(&dA,(size_t)Mc*16)); CK(cudaMalloc(&dC,(size_t)Mc*16));
  fillF<<<(M+TPB-1)/TPB,TPB>>>(dF,N,M,h); CK(cudaDeviceSynchronize());

  cufftHandle planF,planI; FK(cufftPlan1d(&planF,M,CUFFT_D2Z,1)); FK(cufftPlan1d(&planI,M,CUFFT_Z2D,1));

  cudaEvent_t t0,t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
  auto timeDirect=[&](float& ms){ CK(cudaEventRecord(t0));
    directConv<<<(K+TPB-1)/TPB,TPB>>>(dF,dCd,N,K); CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
    CK(cudaEventElapsedTime(&ms,t0,t1)); return 0; };
  auto timeFFT=[&](float& ms){ CK(cudaEventRecord(t0));
    FK(cufftExecD2Z(planF,dF,dA));
    cSquare<<<(Mc+TPB-1)/TPB,TPB>>>(dA,dC,Mc);
    FK(cufftExecZ2D(planI,dC,dOut));
    scaleR<<<(M+TPB-1)/TPB,TPB>>>(dOut,M,1.0/M);     // cuFFT C2R is unnormalized
    CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1)); CK(cudaEventElapsedTime(&ms,t0,t1)); return 0; };

  float dms,fms,tmp; timeDirect(tmp); timeFFT(tmp);   // warm up
  float dbest=1e30f,fbest=1e30f;
  for(int r=0;r<5;++r){ if(timeDirect(dms))return 1; if(dms<dbest)dbest=dms;
                        if(timeFFT(fms))return 1;    if(fms<fbest)fbest=fms; }

  std::vector<double> Cd(K),Cf(M);
  CK(cudaMemcpy(Cd.data(),dCd,(size_t)K*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(Cf.data(),dOut,(size_t)M*8,cudaMemcpyDeviceToHost));
  double cmax=0; for(int k=0;k<K;++k) cmax=fmax(cmax,fabs(Cd[k]));
  double maxRel=0; for(int k=0;k<K;++k) if(fabs(Cd[k])>1e-9*cmax)
    maxRel=fmax(maxRel,fabs(Cf[k]-Cd[k])/fabs(Cd[k]));

  printf("Parton-luminosity convolution L(eta)=int_0^eta F(xi)F(eta-xi)dxi  (N=%d, FFT length M=%d)\n",N,M);
  printf("  FFT vs direct max relerr = %.2e   (same discrete convolution, FFT in O(N log N))\n",maxRel);
  printf("  direct O(N^2) = %.3f ms   |   FFT O(N log N) = %.3f ms   |   speedup = %.1fx\n",
         dbest,fbest,dbest/fbest);
  bool ok=(maxRel<1e-5)&&(fbest<dbest);  // ~1e-7 is the expected FP64 FFT-convolution roundoff floor
  printf("VALIDATION: %s (FFT reproduces the convolution and is faster)\n", ok?"PASS":"FAIL");
  printf("NOTE: applies to the INCLUSIVE parton luminosity only; the event generator must\n"
         "      sample exclusive (x1,x2) per event, which a transform cannot provide.\n");
  cufftDestroy(planF); cufftDestroy(planI);
  cudaFree(dF);cudaFree(dCd);cudaFree(dOut);cudaFree(dA);cudaFree(dC);
  return ok?0:2;
}
