// Validate the Dalitz matrix-element shapes (-DDALITZ_ME) in decay_inc.cuh. Decays many omega(223)
// and eta(221) AT REST, keeps the pi+pi-pi0 final states, and measures the Dalitz-plot density shape:
//   omega/phi P-wave  |M|^2 ~ |p+ x p-|^2  -> vanishes at the Dalitz boundary (edge-suppressed);
//                         report <|p+ x p-|^2>/wmax (higher for ME) and the low-bin (edge) fraction.
//   eta linear slope  |M|^2 ~ 1 + a y + ... (a<0) -> weight pushed to y<0; report <y> (negative for ME).
// Compile WITHOUT the flag -> flat phase space (the control); WITH -DDALITZ_ME -> the ME shape. Run both
// and compare: ME must edge-suppress omega (lower edge fraction, higher mean) and give <y_eta> < 0.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 [-DDALITZ_ME] -DDECAYS -o dalitz_test dalitz_test.cu
// Run:   ./dalitz_test [nPerSpecies=2000000]
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"
#define CK(c) do{cudaError_t e_=(c); if(e_!=cudaSuccess){printf("CUDA %s @ %s:%d\n",cudaGetErrorString(e_),__FILE__,__LINE__);return 1;}}while(0)
#include "decay_inc.cuh"

// Decay one parent (223 or 221) at rest; if the finals are exactly {pi+,pi-,pi0}, return the two Dalitz
// observables: w = |p+ x p-|^2 / wmax  (omega P-wave variable, in [0,1]) and y = 3 T_pi0 / Q - 1.
__host__ __device__ inline bool oneDalitz(uint64_t base,int e,int parent,double& w,double& y){
  // parent generated AT REST (no setup RNG); the decay is driven solely by the dctr stream below.
  double M=dPoleMass(parent);
  double P[4]={0,0,0,M}; int id[1]={parent};
  double F[MAXFINAL*4]; int fid[MAXFINAL];
  uint64_t dctr=(base^0xDA112ULL)+(uint64_t)e*0x100000001B3ULL;
  int nF=decayEvent(P,id,1,dctr,F,fid);
  if(nF!=3) return false;
  int ip=-1,im=-1,i0=-1;
  for(int i=0;i<3;++i){ if(fid[i]==211)ip=i; else if(fid[i]==-211)im=i; else if(fid[i]==111)i0=i; }
  if(ip<0||im<0||i0<0) return false;
  double mpi=0.13957, mpi0=0.13498;
  // |p+ x p-|^2
  double a0=F[4*ip],a1=F[4*ip+1],a2=F[4*ip+2], b0=F[4*im],b1=F[4*im+1],b2=F[4*im+2];
  double cx=a1*b2-a2*b1, cy=a2*b0-a0*b2, cz=a0*b1-a1*b0, cr=cx*cx+cy*cy+cz*cz;
  double M2=M*M, u23=mpi+mpi0, u13=mpi+mpi0;
  double l1=(M2-(mpi+u23)*(mpi+u23))*(M2-(mpi-u23)*(mpi-u23));
  double l2=(M2-(mpi+u13)*(mpi+u13))*(M2-(mpi-u13)*(mpi-u13));
  double wmax=(0.25*fmax(0.0,l1)/M2)*(0.25*fmax(0.0,l2)/M2);
  w=(wmax>0)?cr/wmax:0.0;
  double Q=M-2.0*mpi-mpi0, T0=F[4*i0+3]-mpi0;
  y=3.0*T0/Q-1.0;
  return true;
}

__global__ void kern(int N,uint64_t base,int parent,double* wo,double* yo,int* ok){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  double w,y; bool g=oneDalitz(base,e,parent,w,y); ok[e]=g?1:0; wo[e]=g?w:0; yo[e]=g?y:0;
}

static void stats(const std::vector<double>&w,const std::vector<double>&y,const std::vector<int>&ok,
                  const char* name){
  double sw=0,sy=0; long n=0,edge=0;
  for(size_t i=0;i<ok.size();++i){ if(!ok[i])continue; n++; sw+=w[i]; sy+=y[i]; if(w[i]<0.02) edge++; }
  if(n==0){ printf("  %-6s: no pi+pi-pi0 finals\n",name); return; }
  printf("  %-6s: kept %ld   <|p+xp-|^2/wmax>=%.4f   edge(w<0.02) frac=%.4f   <y>=%+.4f\n",
         name,n,sw/n,(double)edge/n,sy/n);
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):2000000; int TPB=128, blocks=(N+TPB-1)/TPB; uint64_t base=0xDA1A;
  double *dW,*dY; int *dOk;
  CK(cudaMalloc(&dW,(size_t)N*8)); CK(cudaMalloc(&dY,(size_t)N*8)); CK(cudaMalloc(&dOk,(size_t)N*4));
  std::vector<double> hW(N),hY(N); std::vector<int> hOk(N);
#ifdef DALITZ_ME
  printf("Dalitz ME-shape test (-DDALITZ_ME ON), %d decays/species:\n",N);
#else
  printf("Dalitz ME-shape test (FLAT control, no flag), %d decays/species:\n",N);
#endif
  int sp[2]={223,221}; const char* nm[2]={"omega","eta"};
  for(int s=0;s<2;++s){
    kern<<<blocks,TPB>>>(N,base,sp[s],dW,dY,dOk); CK(cudaDeviceSynchronize());
    CK(cudaMemcpy(hW.data(),dW,(size_t)N*8,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(hY.data(),dY,(size_t)N*8,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(hOk.data(),dOk,(size_t)N*4,cudaMemcpyDeviceToHost));
    stats(hW,hY,hOk,nm[s]);
  }
  cudaFree(dW);cudaFree(dY);cudaFree(dOk);
  return 0;
}
