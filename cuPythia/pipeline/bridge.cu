// cuPythia pipeline — shower->hadronization BRIDGE (end-to-end e+e- -> hadrons).
//
// Runs the GPU FSR dipole shower (shower_inc.cuh, the validated core) and writes each event's
// colour-ordered parton chain (q, gluons, qbar = a gluon-KINKED string) to a file, with a
// valid open-singlet colour flow assigned from the chain order. bridge_pythia.cc then reads it
// and hadronizes the kinked string with Pythia's StringFragmentation (forceHadronLevel),
// producing complete hadron-level events. This CLOSES the chain hard->shower->hadronization and
// validates that the GPU shower emits a physically valid, hadronizable colour singlet. (An
// all-GPU gluon-kinked hadronizer — multi-region strings — is the documented next step; the
// GPU hadronize.cu currently does straight q-qbar strings only.)
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o bridge bridge.cu
// Run:   ./bridge [nEvents=5000] [out=shower_partons.dat]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cuda_runtime.h>
#include "shower_inc.cuh"

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)

__global__ void kern(int N,uint64_t base,double* Pout,int* idout,int* nout){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  double P[MAXP*4]; int id[MAXP];
  int n=showerEvent(P,id, base+(uint64_t)e*0x9E3779B97F4A7C15ULL);
  nout[e]=n;
  for(int i=0;i<n;++i){ for(int k=0;k<4;++k) Pout[(size_t)e*MAXP*4+4*i+k]=P[4*i+k]; idout[(size_t)e*MAXP+i]=id[i]; }
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):5000;
  const char* path=(argc>2)?argv[2]:"shower_partons.dat";
  int TPB=128, blocks=(N+TPB-1)/TPB;
  double* dP; int *dId,*dN;
  CK(cudaMalloc(&dP,(size_t)N*MAXP*4*8)); CK(cudaMalloc(&dId,(size_t)N*MAXP*4)); CK(cudaMalloc(&dN,(size_t)N*4));
  kern<<<blocks,TPB>>>(N,0x5110ULL,dP,dId,dN); CK(cudaDeviceSynchronize());
  std::vector<double> hP((size_t)N*MAXP*4); std::vector<int> hId((size_t)N*MAXP),hN(N);
  CK(cudaMemcpy(hP.data(),dP,(size_t)N*MAXP*4*8,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hId.data(),dId,(size_t)N*MAXP*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hN.data(),dN,(size_t)N*4,cudaMemcpyDeviceToHost));

  // Write the chains with an open-singlet colour flow: tag[k]=101+k shared between parton k
  // and k+1, so q(col=101,acol=0), g_k(col=101+k,acol=100+k), qbar(col=0,acol=101+n-2).
  FILE* f=fopen(path,"w"); if(!f){ printf("cannot open %s\n",path); return 1; }
  fprintf(f,"%d %.8f\n",N,MZ);
  for(int e=0;e<N;++e){ int n=hN[e]; fprintf(f,"%d\n",n);
    for(int i=0;i<n;++i){ int id=hId[(size_t)e*MAXP+i];
      int col = (id<0)?0:101+i;            // qbar has no colour; q & gluons do (tag i)
      int acol= (i==0)?0:100+i;            // only the q endpoint has no anticolour; gluons & qbar do (tag i-1)
      const double* p=&hP[(size_t)e*MAXP*4+4*i];
      fprintf(f,"%d %d %d %.8e %.8e %.8e %.8e\n",id,col,acol,p[0],p[1],p[2],p[3]); }
  }
  fclose(f);
  printf("wrote %d GPU-shower parton chains (gluon-kinked singlets) to %s\n",N,path);
  cudaFree(dP);cudaFree(dId);cudaFree(dN);
  return 0;
}
