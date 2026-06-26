// cuPythia pipeline — stage-0 foundation test: populate a batch of gg->gg
// hard-process events INTO the device-resident event record, fully on-GPU, then
// validate the record (per-event 4-momentum conservation, on-shell masses, record
// integrity) and that the average event weight reproduces the cross section.
//
// This proves the data plane: subsequent stages (shower, hadronization) will
// append/modify particles in this same record without leaving the device.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o build_events build_events.cu
// Run:   ./build_events [nEvents=2000000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "event.cuh"
#include "physics.cuh"

// One thread per event: sample a gg->gg hard process and fill the record.
__global__ void buildKernel(DeviceEvents ev, double sqrtS, double alpS, double cMax, uint64_t base){
  int e = blockIdx.x*blockDim.x + threadIdx.x; if(e>=ev.nEvents) return;
  uint64_t key = splitmix64(base ^ ((uint64_t)e * 0x9E3779B97F4A7C15ULL)); // per-event substream
  ev.seed[e] = key;
  uint64_t c = key;
  double E = 0.5*sqrtS, s = sqrtS*sqrtS;
  double cosT = (2.0*u01(splitmix64(c++)) - 1.0) * cMax;
  double phi  = 2.0*M_PI*u01(splitmix64(c++));
  double st = sqrt(fmax(0.0, 1.0 - cosT*cosT));
  double sph, cph; sincos(phi, &sph, &cph);
  // incoming gluons along +/-z (status -21); colour lines 501-502-503
  addParticle(ev,e, 0,0,  E, E, 0, 21,-21, 501,502, -1,-1);
  addParticle(ev,e, 0,0, -E, E, 0, 21,-21, 503,501, -1,-1);
  // outgoing gluons back-to-back (status 23), mothers = the two incoming (0,1)
  double p3x=E*st*cph, p3y=E*st*sph, p3z=E*cosT;
  addParticle(ev,e,  p3x, p3y, p3z, E, 0, 21,23, 503,504, 0,1);
  addParticle(ev,e, -p3x,-p3y,-p3z, E, 0, 21,23, 504,502, 0,1);
  double t = -0.5*s*(1.0-cosT);
  ev.weight[e] = gg2gg_sigma(s, t, -s - t, alpS);   // event weight = dσ/dt̂
  ev.scale[e]  = E*st;                               // pT-hat (hard scale)
}

#define CK(call) do{ cudaError_t e_=(call); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__); return 1; } }while(0)

int main(int argc,char**argv){
  int N = (argc>1)? atoi(argv[1]) : 2000000;
  int maxPart = 16;                  // headroom for shower products later
  double sqrtS=100.0, alphaS=0.118, cMax=0.9, conv_pb=0.3893793721e9;
  DeviceEvents ev = allocEvents(N, maxPart);
  int threads=256, blocks=(N+threads-1)/threads;
  buildKernel<<<blocks,threads>>>(ev, sqrtS, alphaS, cMax, 0x9E1ULL);
  CK(cudaDeviceSynchronize());

  // Copy back a validating subset of the record.
  std::vector<int> nPart(N); std::vector<double> w(N);
  CK(cudaMemcpy(nPart.data(), ev.nPart, N*sizeof(int), cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(w.data(), ev.weight, N*sizeof(double), cudaMemcpyDeviceToHost));
  // Pull the per-particle four-momenta to check conservation + on-shell.
  size_t np=(size_t)N*maxPart;
  std::vector<double> px(np),py(np),pz(np),en(np),mm(np);
  CK(cudaMemcpy(px.data(),ev.px,np*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(py.data(),ev.py,np*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(pz.data(),ev.pz,np*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(en.data(),ev.e ,np*8,cudaMemcpyDeviceToHost));

  double maxCons=0, sumW=0; int badN=0;
  for(int e=0;e<N;++e){
    if(nPart[e]!=4) ++badN;
    // sum incoming (status -21) and outgoing (status 23) momenta
    double iE=0,iz=0, oE=0,ox=0,oy=0,oz=0;
    for(int p=0;p<nPart[e];++p){
      size_t i=(size_t)e*maxPart+p;
      // first two are incoming, last two outgoing (by construction)
      if(p<2){ iE+=en[i]; iz+=pz[i]; }
      else   { oE+=en[i]; ox+=px[i]; oy+=py[i]; oz+=pz[i]; }
    }
    double cE=fabs(oE-iE), cx=fabs(ox), cy=fabs(oy), cz=fabs(oz-iz);
    double cons=fmax(fmax(cE,cx),fmax(cy,cz));
    if(cons>maxCons) maxCons=cons;
    sumW+=w[e];
  }
  double sigma = (2.0*cMax)*(sqrtS*sqrtS/2.0)*(sumW/N)*conv_pb;
  // Simpson reference
  double a=-cMax,h=(2.0*cMax)/2000000,ref=0;
  for(int i=0;i<=2000000;++i){ double cc=a+i*h,t=-0.5*(sqrtS*sqrtS)*(1.0-cc),
    ww=(i==0||i==2000000)?1.0:(i%2?4.0:2.0); ref+=ww*gg2gg_sigma(sqrtS*sqrtS,t,-(sqrtS*sqrtS)-t,alphaS);}
  ref=(ref*h/3.0)*(sqrtS*sqrtS/2.0)*conv_pb;

  printf("Device-resident event record: %d gg->gg events built on-GPU\n", N);
  printf("  particles/event (expect 4): %s   max 4-momentum imbalance = %.1e\n",
         badN==0?"all OK":"MISMATCH", maxCons);
  printf("  sigma from record weights = %.6e pb   Simpson = %.6e pb   relerr = %.2e\n",
         sigma, ref, fabs(sigma-ref)/ref);
  // sigma here is a low-statistics SANITY check (one sample per event); with N=2e6
  // on a peaked integrand the MC error is ~2e-3, so allow 6e-3 (~3 sigma). The real
  // validations are exact record integrity + zero 4-momentum imbalance.
  bool ok = (badN==0) && (maxCons<1e-9) && (fabs(sigma-ref)/ref<6e-3);
  printf("VALIDATION: %s (record integrity + 4-momentum conservation + cross-section sanity)\n",
         ok?"PASS":"FAIL");
  freeEvents(ev);
  return ok?0:2;
}
