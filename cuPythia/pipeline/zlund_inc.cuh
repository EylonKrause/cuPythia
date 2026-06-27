// Shared Lund f(z) sampler — exact port of Pythia 8.317 StringZ::zLundMax/zLund/initFunc.
// Validated in isolation against the analytic f(z) by zlund_test.cu (chi2/ndf ~ 1 in all
// three envelope regimes). Included by hadronize.cu; counter-RNG (../common/rng.cuh).
#pragma once
#include <cmath>
#include <cstdint>
#include "../common/rng.cuh"

__host__ __device__ inline double zLundMax(double a,double b,double c){
  const double AFROMZERO=0.02, AFROMC=0.01;
  bool aIsZero=(a<AFROMZERO), aIsC=(fabs(a-c)<AFROMC);
  double zMax;
  if(aIsZero) zMax=(c>b)? b/c : 1.0;
  else if(aIsC) zMax=b/(b+c);
  else { zMax=0.5*(b+c-sqrt((b-c)*(b-c)+4.0*a*b))/(c-a);
         if(zMax>0.9999 && b>100.0) zMax=fmin(zMax,1.0-a/b); }
  return zMax;
}
__host__ __device__ inline double zLundSample(double a,double b,double c,uint64_t& ctr){
  const double CFROMUNITY=0.01, AFROMZERO=0.02, EXPMAX=50.0;
  bool cIsUnity=(fabs(c-1.0)<CFROMUNITY), aIsZero=(a<AFROMZERO);
  double zMax=zLundMax(a,b,c);
  bool peakedNearZero=(zMax<0.1), peakedNearUnity=(zMax>0.85 && b>1.0);
  double fIntLow=1.,fIntHigh=1.,fInt=2.,zDiv=0.5,zDivC=0.5;
  if(peakedNearZero){
    zDiv=2.75*zMax; fIntLow=zDiv;
    if(cIsUnity) fIntHigh=-zDiv*log(zDiv);
    else { zDivC=pow(zDiv,1.0-c); fIntHigh=zDiv*(1.0-1.0/zDivC)/(c-1.0); }
    fInt=fIntLow+fIntHigh;
  } else if(peakedNearUnity){
    double rcb=sqrt(4.0+(c/b)*(c/b));
    zDiv=rcb-1.0/zMax-(c/b)*log(zMax*0.5*(rcb+c/b));
    if(!aIsZero) zDiv+=(a/b)*log(1.0-zMax);
    zDiv=fmin(zMax,fmax(0.0,zDiv));
    fIntLow=1.0/b; fIntHigh=1.0-zDiv; fInt=fIntLow+fIntHigh;
  }
  double z=0.5; bool accept=false;
  do{
    z=u01(splitmix64(ctr++));
    double fPrel=1.0;
    if(peakedNearZero){
      if(fInt*u01(splitmix64(ctr++))<fIntLow) z=zDiv*z;
      else if(cIsUnity){ z=pow(zDiv,z); fPrel=zDiv/z; }
      else { z=pow(zDivC+(1.0-zDivC)*z,1.0/(1.0-c)); fPrel=pow(zDiv/z,c); }
    } else if(peakedNearUnity){
      if(fInt*u01(splitmix64(ctr++))<fIntLow){ z=zDiv+log(z)/b; fPrel=exp(b*(z-zDiv)); }
      else z=zDiv+(1.0-zDiv)*z;
    }
    if(z>0.0 && z<1.0){
      double fRnd=u01(splitmix64(ctr++));
      double fExp=b*(1.0/zMax-1.0/z)+c*log(zMax/z);
      if(!aIsZero) fExp+=a*log((1.0-z)/(1.0-zMax));
      double fVal=exp(fmax(-EXPMAX,fmin(EXPMAX,fExp)));
      accept=((fVal/fPrel)>fRnd);
    }
  } while(!accept);
  return z;
}
