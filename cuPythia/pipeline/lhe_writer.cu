// cuPythia pipeline — standard I/O: write device-generated gg->gg hard events to a
// spec-valid Les Houches Event File (LHEF 3.0). Events are generated on the GPU (one per
// thread, counter-RNG) in the partonic CM with EXACT 4-momentum conservation and a VALID
// closed colour flow (g1->g2->g3->g4->g1 large-Nc ring); the host writes the LHE. The file
// is meant to be read straight back by Pythia (lhe_validate.cc / Beams:frameType=4), the
// canonical interface that lets cuPythia feed any standard shower/detector toolchain.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o lhe_writer lhe_writer.cu
// Run:   ./lhe_writer [nEvents=10000] [out=events.lhe]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)

// Per event: sample sqrt(sHat) and the CM polar angle (pT-hat cut), build the outgoing
// gluon p3 = (px,py,pz) [p4 = -p3, both energy E]; incoming are (0,0,+-E,E).
__global__ void gen(int N,uint64_t base,double sMin,double sMax,double pTmin,double* out){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  uint64_t ctr=base+(uint64_t)e*0x9E3779B97F4A7C15ULL;
  double shat,E,pT,cth,phi;
  for(int it=0;it<10000;++it){
    shat=sMin+(sMax-sMin)*u01(splitmix64(ctr++)); E=0.5*sqrt(shat);
    cth=2.0*u01(splitmix64(ctr++))-1.0; pT=E*sqrt(fmax(0.0,1.0-cth*cth));
    if(pT>=pTmin) break;
  }
  phi=2.0*M_PI*u01(splitmix64(ctr++));
  out[5*e+0]=E; out[5*e+1]=E*sqrt(fmax(0.0,1.0-cth*cth))*cos(phi);
  out[5*e+2]=E*sqrt(fmax(0.0,1.0-cth*cth))*sin(phi); out[5*e+3]=E*cth; out[5*e+4]=shat;
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):10000;
  const char* path=(argc>2)?argv[2]:"events.lhe";
  double Ebeam=6500.0, sMin=200.0*200.0, sMax=2000.0*2000.0, pTmin=50.0, xsec_pb=1.0e3;
  int TPB=128, blocks=(N+TPB-1)/TPB;
  double* dO; CK(cudaMalloc(&dO,(size_t)N*5*8));
  gen<<<blocks,TPB>>>(N,0x13E5ULL,sMin,sMax,pTmin,dO); CK(cudaDeviceSynchronize());
  std::vector<double> h((size_t)N*5); CK(cudaMemcpy(h.data(),dO,(size_t)N*5*8,cudaMemcpyDeviceToHost));

  FILE* f=fopen(path,"w");
  if(!f){ printf("cannot open %s\n",path); return 1; }
  fprintf(f,"<LesHouchesEvents version=\"3.0\">\n");
  fprintf(f,"<header>\n  cuPythia GPU gg->gg parton-level events (partonic CM, pT-hat>%.0f GeV)\n</header>\n",pTmin);
  fprintf(f,"<init>\n");
  // IDBMUP1 IDBMUP2 EBMUP1 EBMUP2 PDFGUP1 PDFGUP2 PDFSUP1 PDFSUP2 IDWTUP NPRUP
  fprintf(f,"2212 2212 %.8E %.8E 0 0 0 0 3 1\n",Ebeam,Ebeam);
  // XSECUP XERRUP XMAXUP LPRUP
  fprintf(f,"%.8E %.8E 1.0E+00 1\n",xsec_pb,xsec_pb*0.01);
  fprintf(f,"</init>\n");
  for(int e=0;e<N;++e){
    double E=h[5*e],p3x=h[5*e+1],p3y=h[5*e+2],p3z=h[5*e+3],scale=sqrt(h[5*e+4]);
    fprintf(f,"<event>\n");
    // NUP IDPRUP XWGTUP SCALUP AQEDUP AQCDUP
    fprintf(f,"4 1 1.0E+00 %.8E 7.546771E-03 1.180000E-01\n",scale);
    // id status mo1 mo2 col acol px py pz E m vtim spin  (valid gg->gg t-channel flow,
    // the convention Pythia itself emits: incoming col/acol act swapped in line-tracing)
    fprintf(f,"21 -1 0 0 101 102 0.0 0.0 %.8E %.8E 0.0 0. 1.\n", E, E);
    fprintf(f,"21 -1 0 0 103 101 0.0 0.0 %.8E %.8E 0.0 0. 1.\n",-E, E);
    fprintf(f,"21  1 1 2 103 104 %.8E %.8E %.8E %.8E 0.0 0. 1.\n", p3x, p3y, p3z, E);
    fprintf(f,"21  1 1 2 104 102 %.8E %.8E %.8E %.8E 0.0 0. 1.\n",-p3x,-p3y,-p3z, E);
    fprintf(f,"</event>\n");
  }
  fprintf(f,"</LesHouchesEvents>\n");
  fclose(f);
  printf("wrote %d gg->gg events to %s (spec-valid LHEF 3.0, valid t-channel colour flow)\n",N,path);
  cudaFree(dO);
  return 0;
}
