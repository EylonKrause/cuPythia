// cuPythia kernel 12 — 2->2 phase-space generation (full final-state kinematics).
//
// Kernels 02/03 sampled only the scattering angle. Here each GPU thread generates
// the full final-state FOUR-MOMENTA of a 2->2 event (what you would write to an
// LHE / HepMC record), for both a massless (gg->gg) and a massive (m=1.5 GeV,
// charm-like) final state. Validated by:
//   * exact energy-momentum conservation   p3 + p4 = p1 + p2   (machine precision)
//   * on-shell masses   p3^2 = p4^2 = m^2
//   * invariant mass    (p3 + p4)^2 = s
//   * massless cross section reproduced vs Simpson quadrature.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o phase_space phase_space.cu
// Run:   ./phase_space [trialsPerThread=4000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

struct V4 { double e, x, y, z; };
__host__ __device__ inline double dot4(const V4& a, const V4& b) {
  return a.e*b.e - a.x*b.x - a.y*b.y - a.z*b.z;
}
__host__ __device__ inline double pow2(double x){ return x*x; }
__host__ __device__ inline double gg2gg_sigma(double s,double t,double u,double aS){
  double is=1.0/s,it=1.0/t,iu=1.0/u; // reciprocal precompute: 3 FP64 div instead of 13
  double rts=t*is,rst=s*it,rus=u*is,rsu=s*iu,rtu=t*iu,rut=u*it;
  double a=(9./4.)*(rts*rts+2.*rts+3.+2.*rst+rst*rst);
  double b=(9./4.)*(rus*rus+2.*rus+3.+2.*rsu+rsu*rsu);
  double c=(9./4.)*(rtu*rtu+2.*rtu+3.+2.*rut+rut*rut);
  return (M_PI*is*is)*pow2(aS)*0.5*(a+b+c);
}

// Generate 2->2 four-momenta in the CM frame; reduce conservation / mass residuals
// and (for the massless run) the cross-section integrand.
__global__ void psKernel(uint64_t seed, uint64_t nPer, double sqrtS, double m,
                         double cMax, double alpS, int massless,
                         double* gMaxCons, double* gMaxMass, double* gMaxShat, double* gSumSig){
  uint64_t tid = blockIdx.x*(uint64_t)blockDim.x + threadIdx.x;
  uint64_t ctr = seed + tid*0x100000001B3ULL;
  double E = 0.5*sqrtS, s = sqrtS*sqrtS;
  double p = sqrt(fmax(0.0, E*E - m*m));    // final-state momentum magnitude
  V4 p1 = {E, 0, 0, E}, p2 = {E, 0, 0, -E}; // incoming along z
  double lCons=0, lMass=0, lShat=0, lSig=0;
  for (uint64_t i=0;i<nPer;++i){
    double c = (2.0*u01(splitmix64(ctr++))-1.0)*cMax;        // cos theta
    double phi = 2.0*M_PI*u01(splitmix64(ctr++));
    double st = sqrt(fmax(0.0,1.0-c*c));
    V4 p3 = {E,  p*st*cos(phi),  p*st*sin(phi),  p*c};
    V4 p4 = {E, -p*st*cos(phi), -p*st*sin(phi), -p*c};
    // conservation: p3+p4 must equal p1+p2 = (sqrtS,0,0,0)
    double dE=(p3.e+p4.e)-(p1.e+p2.e), dx=p3.x+p4.x, dy=p3.y+p4.y, dz=p3.z+p4.z;
    double cons = fmax(fabs(dE), fmax(fabs(dx), fmax(fabs(dy), fabs(dz))));
    if (cons>lCons) lCons=cons;
    // on-shell: p3^2 = m^2
    double mres = fabs(dot4(p3,p3) - m*m);
    if (mres>lMass) lMass=mres;
    // invariant mass: (p3+p4)^2 = s
    V4 sum={p3.e+p4.e,p3.x+p4.x,p3.y+p4.y,p3.z+p4.z};
    double sres = fabs(dot4(sum,sum) - s);
    if (sres>lShat) lShat=sres;
    // massless cross section from the reconstructed kinematics
    if (massless){
      double t = -0.5*s*(1.0 - c);
      lSig += gg2gg_sigma(s, t, -s - t, alpS);
    }
  }
  atomicMax((unsigned long long*)gMaxCons, (unsigned long long)__double_as_longlong(lCons));
  atomicMax((unsigned long long*)gMaxMass, (unsigned long long)__double_as_longlong(lMass));
  atomicMax((unsigned long long*)gMaxShat, (unsigned long long)__double_as_longlong(lShat));
  atomicAdd(gSumSig, lSig);
}
static double simpson(int N,double s,double aS,double cMax){
  double a=-cMax,h=(2.0*cMax)/N,sum=0;
  for(int i=0;i<=N;++i){ double c=a+i*h,t=-0.5*s*(1.0-c),w=(i==0||i==N)?1.0:(i%2?4.0:2.0);
    sum+=w*gg2gg_sigma(s,t,-s-t,aS);} return (sum*h/3.0)*(s/2.0);
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

static int run(uint64_t seed,uint64_t nPer,int blocks,int threads,double sqrtS,double m,
               double cMax,double alpS,int massless,double* mc,double* mm,double* ms,double* sig){
  double *dC,*dM,*dS,*dSig;
  CK(cudaMalloc(&dC,8)); CK(cudaMalloc(&dM,8)); CK(cudaMalloc(&dS,8)); CK(cudaMalloc(&dSig,8));
  CK(cudaMemset(dC,0,8)); CK(cudaMemset(dM,0,8)); CK(cudaMemset(dS,0,8)); CK(cudaMemset(dSig,0,8));
  psKernel<<<blocks,threads>>>(seed,nPer,sqrtS,m,cMax,alpS,massless,dC,dM,dS,dSig);
  CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(mc,dC,8,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(mm,dM,8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(ms,dS,8,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(sig,dSig,8,cudaMemcpyDeviceToHost));
  cudaFree(dC);cudaFree(dM);cudaFree(dS);cudaFree(dSig); return 0;
}

int main(int argc,char**argv){
  uint64_t nPer=(argc>1)?strtoull(argv[1],nullptr,10):4000ULL;
  const int blocks=1024,threads=256;
  double sqrtS=100.0, alphaS=0.118, cMax=0.9, conv_pb=0.3893793721e9;
  uint64_t total=(uint64_t)blocks*threads*nPer;

  double mc0,mm0,ms0,sig0;  // massless gg->gg
  if(run(0x9501ULL,nPer,blocks,threads,sqrtS,0.0,cMax,alphaS,1,&mc0,&mm0,&ms0,&sig0)) return 1;
  double mc1,mm1,ms1,sig1;  // massive (charm-like m=1.5 GeV)
  if(run(0x9502ULL,nPer,blocks,threads,sqrtS,1.5,cMax,alphaS,0,&mc1,&mm1,&ms1,&sig1)) return 1;

  double sigma = (2.0*cMax)*(sqrtS*sqrtS/2.0)*(sig0/(double)total)*conv_pb;
  double ref = simpson(2000000,sqrtS*sqrtS,alphaS,cMax)*conv_pb;
  printf("2->2 phase-space generation (full four-momenta), sqrt(s)=%.0f GeV\n", sqrtS);
  printf("  massless gg->gg : max|p_cons|=%.1e  max|m^2|=%.1e  max|shat-s|=%.1e\n", mc0, mm0, ms0);
  printf("  massive  (m=1.5): max|p_cons|=%.1e  max|m^2-2.25|=%.1e  max|shat-s|=%.1e\n", mc1, mm1, ms1);
  printf("  cross section (massless) = %.6e pb   Simpson = %.6e pb   relerr=%.2e\n",
         sigma, ref, fabs(sigma-ref)/ref);
  bool ok = mc0<1e-9 && ms0<1e-6 && mc1<1e-9 && mm1<1e-6 && ms1<1e-6 && fabs(sigma-ref)/ref<2e-3;
  printf("VALIDATION: %s (4-momentum conservation + on-shell masses + cross section)\n", ok?"PASS":"FAIL");
  return ok?0:2;
}
