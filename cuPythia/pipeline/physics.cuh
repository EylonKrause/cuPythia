// cuPythia pipeline — shared physics (matrix elements). Reciprocal-optimized.
#pragma once
#include "../common/rng.cuh"   // M_PI fallback

__host__ __device__ inline double pow2(double x){ return x*x; }

// gg -> gg tree-level dσ/dt̂ (verbatim Pythia 8.317 Sigma2gg2gg::sigmaKin),
// reciprocal-precompute: 3 FP64 divisions instead of 13 (~2.9x on Blackwell).
__host__ __device__ inline double gg2gg_sigma(double sH,double tH,double uH,double aS){
  double is=1.0/sH, it=1.0/tH, iu=1.0/uH;
  double rts=tH*is, rst=sH*it, rus=uH*is, rsu=sH*iu, rtu=tH*iu, rut=uH*it;
  double a=(9./4.)*(rts*rts+2.*rts+3.+2.*rst+rst*rst);
  double b=(9./4.)*(rus*rus+2.*rus+3.+2.*rsu+rsu*rsu);
  double c=(9./4.)*(rtu*rtu+2.*rtu+3.+2.*rut+rut*rut);
  return (M_PI*is*is)*pow2(aS)*0.5*(a+b+c);
}
