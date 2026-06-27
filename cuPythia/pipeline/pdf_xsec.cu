// cuPythia pipeline stage 1 — PDF convolution -> a real HADRONIC cross section.
//
// Turns the partonic gg->gg ME into pp->gg+X by convolving with gluon PDFs:
//   sigma = int dx1 dx2 dcos  f_g(x1,muF) f_g(x2,muF) (dsigma_hat/dt_hat)(sHat) (sHat/2),
// with sHat = x1 x2 S and a pT-hat cut. The device PDF (pdf.cuh) is read by
// bilinear interpolation with low-x/Q^2 freezing. Validated two ways: (1) device
// interpolation vs the analytic toy PDF; (2) the device-MC hadronic sigma vs a CPU
// reference on the same RNG samples (convolution machinery correctness).
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o pdf_xsec pdf_xsec.cu
// Run:   ./pdf_xsec [trialsPerThread=4000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "pdf.cuh"
#include "physics.cuh"

__host__ __device__ inline double alphaS1(double mu2,double muRef2,double aSref){
  double b0=(33.0-2.0*5.0)/(12.0*M_PI); return aSref/(1.0 + aSref*b0*log(mu2/muRef2));
}

// One thread does kPer hadronic trials; integrand uses the device PDF grid.
__global__ void hadKernel(PdfGrid g,uint64_t seed,uint64_t kPer,double S,double xmin,
                          double cMax,double sHatMin,double aSref,double muRef2,double* gSum){
  uint64_t tid=blockIdx.x*(uint64_t)blockDim.x+threadIdx.x, ctr=seed+tid*0x100000001B3ULL;
  double local=0.0;
  for(uint64_t i=0;i<kPer;++i){
    double x1=xmin+(1.0-xmin)*u01(splitmix64(ctr++));
    double x2=xmin+(1.0-xmin)*u01(splitmix64(ctr++));
    double cosT=(2.0*u01(splitmix64(ctr++))-1.0)*cMax;
    double sHat=x1*x2*S;
    if(sHat<sHatMin) continue;
    double pT=0.5*sqrt(sHat)*sqrt(fmax(0.0,1.0-cosT*cosT)); double mu0=(pT>1.0)?pT:1.0;
    double aS=alphaS1(mu0*mu0,muRef2,aSref);
    double t=-0.5*sHat*(1.0-cosT);
    double dsig=gg2gg_sigma(sHat,t,-sHat-t,aS);
    double f1=pdf_g(g,x1,mu0*mu0), f2=pdf_g(g,x2,mu0*mu0);
    local += f1*f2*dsig*(sHat/2.0);           // dt_hat = sHat/2 dcos
  }
  atomicAdd(gSum,local);
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  uint64_t kPer=(argc>1)?strtoull(argv[1],nullptr,10):4000ULL;
  int blocks=16,threads=256;                      // 4096 threads: CPU reference mirrors the SAME samples
  double rootS=13000.0, S=rootS*rootS;            // 13 TeV pp
  double xmin=1e-3, cMax=0.8, pTmin=50.0;         // 50 GeV jet pT-hat cut
  double sHatMin=(2.0*pTmin)*(2.0*pTmin)/(1.0-cMax*cMax); // ensure pT-hat>=pTmin reachable
  double aSref=0.118, muRef2=91.1876*91.1876, conv_pb=0.3893793721e9;
  uint64_t seed=0x9D7ULL;

  PdfGrid g = buildToyGrid(200,60,1e-4,1.0,1.0,2000.0);
  std::vector<double> hxf((size_t)g.nx*g.nq); CK(cudaMemcpy(hxf.data(),g.xf,(size_t)g.nx*g.nq*8,cudaMemcpyDeviceToHost));
  PdfGrid cg=g; cg.xf=hxf.data();                 // host-pointer copy for host interpolation

  // (1) pointwise interpolation accuracy, measured in the x-range that carries the
  //     cross section at this pT-cut (x ~ 1e-3..0.3); the (1-x)^5 corner at x->1 is
  //     under-resolved by a 200-pt log grid but contributes ~0 to sigma (see (3)).
  double maxRel=0; uint64_t rc=12345;
  for(int k=0;k<20000;++k){
    double x =exp(log(1e-3)+(log(0.3)-log(1e-3))*u01(splitmix64(rc++)));
    double Q2=exp(log(4.0) +(log(1e6)-log(4.0)) *u01(splitmix64(rc++)));
    double interp=pdf_xfg(cg,x,Q2), exact=toy_xfg(x,Q2);
    maxRel=fmax(maxRel,fabs(interp-exact)/fabs(exact)); }

  // (2) hadronic sigma on the GPU (grid interpolation).
  double* dS; CK(cudaMalloc(&dS,8)); CK(cudaMemset(dS,0,8));
  hadKernel<<<blocks,threads>>>(g,seed,kPer,S,xmin,cMax,sHatMin,aSref,muRef2,dS);
  CK(cudaDeviceSynchronize());
  double hsum=0; CK(cudaMemcpy(&hsum,dS,8,cudaMemcpyDeviceToHost));
  uint64_t total=(uint64_t)blocks*threads*kPer;
  double V=(1.0-xmin)*(1.0-xmin)*(2.0*cMax);
  double sigma_dev=V*(hsum/total)*conv_pb;

  // (3) CPU reference over the IDENTICAL samples, two ways:
  //     csum_i = grid interpolation (machinery/determinism check vs the GPU)
  //     csum_e = analytic toy PDF   (interpolation FIDELITY at the cross-section level)
  double csum_i=0, csum_e=0;
  for(uint64_t tid=0; tid<(uint64_t)blocks*threads; ++tid){ uint64_t ctr=seed+tid*0x100000001B3ULL;
    for(uint64_t i=0;i<kPer;++i){
      double x1=xmin+(1.0-xmin)*u01(splitmix64(ctr++)); double x2=xmin+(1.0-xmin)*u01(splitmix64(ctr++));
      double cosT=(2.0*u01(splitmix64(ctr++))-1.0)*cMax; double sHat=x1*x2*S; if(sHat<sHatMin) continue;
      double pT=0.5*sqrt(sHat)*sqrt(fmax(0.0,1.0-cosT*cosT)); double mu0=(pT>1.0)?pT:1.0;
      double aS=alphaS1(mu0*mu0,muRef2,aSref); double t=-0.5*sHat*(1.0-cosT);
      double common=gg2gg_sigma(sHat,t,-sHat-t,aS)*(sHat/2.0);
      csum_i += pdf_g(cg,x1,mu0*mu0)*pdf_g(cg,x2,mu0*mu0)*common;                  // grid interp
      csum_e += (toy_xfg(x1,mu0*mu0)/x1)*(toy_xfg(x2,mu0*mu0)/x2)*common; } }       // analytic
  double sigma_cpu_i=V*(csum_i/total)*conv_pb;     // same sample count as the GPU
  double sigma_cpu_e=V*(csum_e/total)*conv_pb;

  double rel_machinery=fabs(sigma_dev-sigma_cpu_i)/sigma_cpu_i;  // GPU vs CPU, same samples
  double rel_interp   =fabs(sigma_cpu_i-sigma_cpu_e)/sigma_cpu_e;// grid vs analytic, same samples

  printf("PDF convolution -> hadronic gg->gg cross section (13 TeV pp, pT-hat>%.0f GeV)\n",pTmin);
  printf("  pointwise interp relerr in sigma-support (x:1e-3..0.3) = %.2e\n", maxRel);
  printf("  hadronic sigma (GPU, grid interp)      = %.4e pb\n", sigma_dev);
  printf("  hadronic sigma (CPU, grid interp)      = %.4e pb   GPU-vs-CPU relerr = %.2e (determinism)\n",
         sigma_cpu_i, rel_machinery);
  printf("  hadronic sigma (CPU, analytic PDF)     = %.4e pb   interp fidelity  = %.2e\n",
         sigma_cpu_e, rel_interp);
  bool ok = (maxRel<3e-2) && (rel_machinery<1e-3) && (rel_interp<1e-2);
  printf("VALIDATION: %s (interp fidelity at sigma level + GPU/CPU determinism)\n", ok?"PASS":"FAIL");
  freeGrid(g); cudaFree(dS);
  return ok?0:2;
}
