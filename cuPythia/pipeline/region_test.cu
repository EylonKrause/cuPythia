// cuPythia — validate the multi-region string-table construction (region_inc.cuh) in
// isolation, before wiring kinematicsHadron into the fragmentation chain (the plan's
// "unit-test region construction standalone first"). Checks, for a gluon-kinked q-g-qbar
// and q-g-g-qbar string: every non-empty region has lightlike pPos/pNeg, an orthonormal
// transverse basis (eX^2=eY^2=-1, eX.pPos=eX.pNeg=eY.pPos=eY.pNeg=eX.eY=0), w2=2 pPos.pNeg,
// and project/pHad are exact inverses (the completeness check). Plus host==device.
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o region_test region_test.cu
// Run:   ./region_test

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <cuda_runtime.h>
#include "../common/rng.cuh"
#include "region_inc.cuh"

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)

// max violation of region math properties + projection round-trip, scaled by sqrt(w2).
__host__ __device__ inline double checkRegion(const Region& R,uint64_t& ctr){
  if(!R.isSetUp || R.isEmpty) return 0.0;
  double s=sqrt(R.w2), v=0.0;
  v=fmax(v,fabs(v4dot(R.pPos,R.pPos))/R.w2);      // lightlike
  v=fmax(v,fabs(v4dot(R.pNeg,R.pNeg))/R.w2);
  v=fmax(v,fabs(v4dot(R.eX,R.eX)+1.0));            // eX^2 = -1
  v=fmax(v,fabs(v4dot(R.eY,R.eY)+1.0));
  v=fmax(v,fabs(v4dot(R.eX,R.pPos))/s);            // orthogonality
  v=fmax(v,fabs(v4dot(R.eX,R.pNeg))/s);
  v=fmax(v,fabs(v4dot(R.eY,R.pPos))/s);
  v=fmax(v,fabs(v4dot(R.eY,R.pNeg))/s);
  v=fmax(v,fabs(v4dot(R.eX,R.eY)));
  v=fmax(v,fabs(R.w2-2.0*v4dot(R.pPos,R.pNeg))/R.w2);
  // project(pHad(a,b,c,d)) == (a,b,c,d) for random coords, and pHad(project(p))==p.
  for(int t=0;t<8;++t){
    double a=2.0*u01(splitmix64(ctr++))-1.0, b=2.0*u01(splitmix64(ctr++))-1.0;
    double c=s*(2.0*u01(splitmix64(ctr++))-1.0), d=s*(2.0*u01(splitmix64(ctr++))-1.0);
    double q[4]; regionPHad(R,a,b,c,d,q);
    double xa,xb,xc,xd; regionProject(R,q,xa,xb,xc,xd);
    v=fmax(v,fabs(xa-a)); v=fmax(v,fabs(xb-b)); v=fmax(v,(fabs(xc-c)+fabs(xd-d))/s);
  }
  return v;
}

// Build a gluon-kinked chain's regions (incl. lazy cross regions for n<=4) and check all.
__host__ __device__ inline double checkChain(const double* P,const int* id,int n,uint64_t ctr){
  StringSys S; sysSetUp(S,P,id,n);
  // lazily build cross regions from regionLowPos/regionLowNeg (isMassless=true)
  for(int iPos=0;iPos<S.sizeStr;++iPos) for(int iNeg=0;iNeg<S.sizeStr;++iNeg){
    if(iPos+iNeg>=S.sizeStr) continue;            // low/valid regions only
    int r=sysIReg(S,iPos,iNeg);
    if(!S.reg[r].isSetUp){
      const Region& lp=S.reg[sysIReg(S,iPos,S.iMax-iPos)];   // regionLowPos
      const Region& ln=S.reg[sysIReg(S,S.iMax-iNeg,iNeg)];   // regionLowNeg
      regionSetUp(S.reg[r], lp.pPos, ln.pNeg, 101,101, true);
    }
  }
  double v=0.0; for(int r=0;r<S.nReg;++r) v=fmax(v,checkRegion(S.reg[r],ctr));
  return v;
}

__global__ void kern(const double* P,const int* id,int n,double* out){
  *out=checkChain(P,id,n,0xC0FFEEULL);
}

int main(){
  double Q=91.1876;
  // q-g-qbar Mercedes (3 massless, sum=(0,0,0,Q)); angles 0/120/240 in x-z.
  double E3=Q/3.0, c=0.8660254037844387;
  double P3[12]={ 0,0,E3,E3,  E3*c,0,-0.5*E3,E3,  -E3*c,0,-0.5*E3,E3 };
  int    id3[3]={1,21,-1};
  // q-g-g-qbar (4 massless, sum=(0,0,0,Q)); angles 45/135/225/315 in x-z.
  double E4=Q/4.0, r=0.7071067811865476;
  double P4[16]={ E4*r,0,E4*r,E4,  E4*r,0,-E4*r,E4,  -E4*r,0,-E4*r,E4,  -E4*r,0,E4*r,E4 };
  int    id4[4]={1,21,21,-1};

  double v3=checkChain(P3,id3,3,0x111ULL);
  double v4=checkChain(P4,id4,4,0x222ULL);

  // device determinism: build+check on GPU, compare to host value for n=3.
  double *dP,*dOut; int* dId;
  CK(cudaMalloc(&dP,12*8)); CK(cudaMalloc(&dId,3*4)); CK(cudaMalloc(&dOut,8));
  CK(cudaMemcpy(dP,P3,12*8,cudaMemcpyHostToDevice)); CK(cudaMemcpy(dId,id3,3*4,cudaMemcpyHostToDevice));
  double hostV3=checkChain(P3,id3,3,0xC0FFEEULL);
  kern<<<1,1>>>(dP,dId,3,dOut); CK(cudaDeviceSynchronize());
  double devV3; CK(cudaMemcpy(&devV3,dOut,8,cudaMemcpyDeviceToHost));

  printf("Multi-region string-table construction (gluon-kinked) — math validation\n");
  printf("  q-g-qbar   (3 regions): max property/round-trip violation = %.2e\n", v3);
  printf("  q-g-g-qbar (6 regions): max property/round-trip violation = %.2e\n", v4);
  printf("  host vs device (n=3)  : |diff| = %.2e\n", fabs(hostV3-devV3));
  bool ok=(v3<1e-9)&&(v4<1e-9)&&(fabs(hostV3-devV3)<1e-10);
  printf("VALIDATION: %s (lightlike basis + orthonormal axes + project/pHad inverse + host==device)\n",
         ok?"PASS":"FAIL");
  cudaFree(dP);cudaFree(dId);cudaFree(dOut);
  return ok?0:2;
}
