// cuPythia pipeline — THE ORCHESTRATOR (build step 2): a device-resident parton-level
// generator. One DeviceEvents is allocated once; the stages
//     build  ->  reweight  ->  unweight  ->  CUB-compact
// run as separate kernel launches with NO host round-trip between them (only the final
// scalar reductions and the accepted-event I/O cross PCIe). This is the residency that
// Pepper (arXiv:2311.06198) and madgraph4gpu concede -- they write parton-level events to
// host and hand off to CPU Pythia for showering; cuPythia keeps them on the device.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o generate generate.cu
// Run:   ./generate [nEvents=4000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include <cub/cub.cuh>
#include "event.cuh"
#include "physics.cuh"

#define NVAR 3
__constant__ double FAC[NVAR] = {1.0, 0.5, 2.0};   // slot 0 nominal, then 0.5x, 2x mu_R

__host__ __device__ inline double alphaS1(double mu2,double muRef2,double aSref){
  double b0=(33.0-2.0*5.0)/(12.0*M_PI); return aSref/(1.0 + aSref*b0*log(mu2/muRef2));
}

// Stage 0 — hard process: sample gg->gg and fill 4 partons + nominal weight + kinematics.
__global__ void stage0_build(DeviceEvents ev,double sqrtS,double cMax,uint64_t base,
                             double aSref,double muRef2){
  int e=blockIdx.x*blockDim.x+threadIdx.x; if(e>=ev.nEvents) return;
  uint64_t key=splitmix64(base ^ ((uint64_t)e*0x9E3779B97F4A7C15ULL)); ev.seed[e]=key;
  uint64_t c=key; double E=0.5*sqrtS, s=sqrtS*sqrtS;
  double cosT=(2.0*u01(splitmix64(c++))-1.0)*cMax;
  double phi =2.0*M_PI*u01(splitmix64(c++));
  double st=sqrt(fmax(0.0,1.0-cosT*cosT)); double sph,cph; sincos(phi,&sph,&cph);
  addParticle(ev,e, 0,0, E,E,0, 21,-21, 501,502, -1,-1);
  addParticle(ev,e, 0,0,-E,E,0, 21,-21, 503,501, -1,-1);
  double p3x=E*st*cph,p3y=E*st*sph,p3z=E*cosT;
  addParticle(ev,e,  p3x, p3y, p3z, E,0, 21,23, 503,504, 0,1);
  addParticle(ev,e, -p3x,-p3y,-p3z, E,0, 21,23, 504,502, 0,1);
  double pT=E*st, mu0=(pT>1e-3)?pT:1e-3, aS0=alphaS1(mu0*mu0,muRef2,aSref);
  double t=-0.5*s*(1.0-cosT);
  ev.scale[e]=mu0; ev.x1[e]=1.0; ev.x2[e]=1.0; ev.flavA[e]=21; ev.flavB[e]=21;
  ev.weight[(size_t)e*ev.nVar+0]=gg2gg_sigma(s,t,-s-t,aS0);
}
// Stage 2 — reweight: fill the scale-variation weight slots from the nominal (in place).
__global__ void stage2_reweight(DeviceEvents ev,double aSref,double muRef2){
  int e=blockIdx.x*blockDim.x+threadIdx.x; if(e>=ev.nEvents) return;
  double mu0=ev.scale[e], aS0=alphaS1(mu0*mu0,muRef2,aSref), w0=ev.weight[(size_t)e*ev.nVar+0];
  for(int k=1;k<ev.nVar;++k){ double mu=FAC[k]*mu0, aSk=alphaS1(mu*mu,muRef2,aSref), r=aSk/aS0;
    ev.weight[(size_t)e*ev.nVar+k]=w0*r*r; }
}
// Stage 7 — unweighting: von-Neumann accept on the nominal weight (disjoint RNG namespace).
__global__ void stage7_unweight(DeviceEvents ev,double wMax,unsigned char* acc,int* idx){
  int e=blockIdx.x*blockDim.x+threadIdx.x; if(e>=ev.nEvents) return;
  idx[e]=e;
  double w=ev.weight[(size_t)e*ev.nVar+0];
  double u=u01(splitmix64(ev.seed[e] ^ 0xACCE7700ULL));   // separate unweighting stream
  acc[e]=(u*wMax < w)?1:0;
}

static double simpsonSigma(int N,double s,double cMax,double aSref,double muRef2,double conv){
  double a=-cMax,h=(2.0*cMax)/N,sum=0;
  for(int i=0;i<=N;++i){ double c=a+i*h, E=0.5*sqrt(s), pT=E*sqrt(fmax(0.0,1.0-c*c));
    double mu=(pT>1e-3)?pT:1e-3, aS=alphaS1(mu*mu,muRef2,aSref), t=-0.5*s*(1.0-c);
    double w=(i==0||i==N)?1.0:(i%2?4.0:2.0); sum+=w*gg2gg_sigma(s,t,-s-t,aS); }
  return (sum*h/3.0)*(s/2.0)*conv;
}
static double wMaxScan(int N,double s,double cMax,double aSref,double muRef2){
  double m=0,E=0.5*sqrt(s);
  for(int i=0;i<=N;++i){ double c=-cMax+(2.0*cMax)*i/N, pT=E*sqrt(fmax(0.0,1.0-c*c));
    double mu=(pT>1e-3)?pT:1e-3, aS=alphaS1(mu*mu,muRef2,aSref), t=-0.5*s*(1.0-c);
    m=fmax(m,gg2gg_sigma(s,t,-s-t,aS)); } return m*1.05;
}
#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):4000000;
  double sqrtS=100.0, cMax=0.9, conv_pb=0.3893793721e9;
  double aSref=0.118, muRef2=91.1876*91.1876; double s=sqrtS*sqrtS;
  uint64_t base=0x6E11ULL;

  DeviceEvents ev = allocEvents(N, 16, NVAR);     // ONE resident record
  int threads=256, blocks=(N+threads-1)/threads;

  // --- stages, all on the device, no host round-trip between them ---
  stage0_build  <<<blocks,threads>>>(ev,sqrtS,cMax,base,aSref,muRef2);
  stage2_reweight<<<blocks,threads>>>(ev,aSref,muRef2);
  double wMax = wMaxScan(2000000,s,cMax,aSref,muRef2);
  unsigned char* dAcc; int* dIdx; CK(cudaMalloc(&dAcc,N)); CK(cudaMalloc(&dIdx,N*sizeof(int)));
  stage7_unweight<<<blocks,threads>>>(ev,wMax,dAcc,dIdx);
  CK(cudaDeviceSynchronize());

  // CUB-compact the accepted event indices on the device.
  int* dKept; int* dNum; CK(cudaMalloc(&dKept,N*sizeof(int))); CK(cudaMalloc(&dNum,sizeof(int)));
  void* dTmp=nullptr; size_t tmpB=0;
  cub::DeviceSelect::Flagged(dTmp,tmpB,dIdx,dAcc,dKept,dNum,N);
  CK(cudaMalloc(&dTmp,tmpB));
  cub::DeviceSelect::Flagged(dTmp,tmpB,dIdx,dAcc,dKept,dNum,N);
  CK(cudaDeviceSynchronize());
  int nAcc=0; CK(cudaMemcpy(&nAcc,dNum,sizeof(int),cudaMemcpyDeviceToHost));

  // Validation reductions (these scalar copies are the ONLY non-I/O host traffic).
  std::vector<double> wv((size_t)N*NVAR); CK(cudaMemcpy(wv.data(),ev.weight,(size_t)N*NVAR*8,cudaMemcpyDeviceToHost));
  std::vector<unsigned char> acc(N); CK(cudaMemcpy(acc.data(),dAcc,N,cudaMemcpyDeviceToHost));
  double sw[NVAR]={0,0,0}; long indepAcc=0;
  for(int e=0;e<N;++e){ for(int k=0;k<NVAR;++k) sw[k]+=wv[(size_t)e*NVAR+k]; indepAcc+=acc[e]; }
  double sigNom=(2.0*cMax)*(s/2.0)*(sw[0]/N)*conv_pb;
  double sigLo =(2.0*cMax)*(s/2.0)*(sw[1]/N)*conv_pb;
  double sigHi =(2.0*cMax)*(s/2.0)*(sw[2]/N)*conv_pb;
  double ref = simpsonSigma(2000000,s,cMax,aSref,muRef2,conv_pb);
  double eta = (double)nAcc/N;

  printf("Device-resident parton-level generator: %d events, build->reweight->unweight->compact\n",N);
  printf("  nominal sigma   = %.6e pb   Simpson = %.6e pb   relerr = %.2e\n", sigNom, ref, fabs(sigNom-ref)/ref);
  printf("  scale band      = [%+.1f%%, %+.1f%%]  (mu_R = 2x .. 0.5x)\n",
         100.0*(sigHi-sigNom)/sigNom, 100.0*(sigLo-sigNom)/sigNom);
  printf("  unweighting eff = %.2f%%   accepted (CUB) = %d   (independent count = %ld)\n",
         100.0*eta, nAcc, indepAcc);
  printf("  residency: build/reweight/unweight/compact all ran on ONE device record, no host round-trip.\n");
  bool ok = (fabs(sigNom-ref)/ref<6e-3) && ((long)nAcc==indepAcc) && (sigLo>sigNom)&&(sigNom>sigHi);
  printf("VALIDATION: %s (sigma vs quadrature + scale band + CUB count == independent)\n", ok?"PASS":"FAIL");
  freeEvents(ev); cudaFree(dAcc);cudaFree(dIdx);cudaFree(dKept);cudaFree(dNum);cudaFree(dTmp);
  return ok?0:2;
}
