// cuPythia hadronization — Breit-Wigner resonance masses for vector mesons (precision
// correction toward theory; pole masses are the default). Samples the vector mass from a
// relativistic Breit-Wigner (Lorentzian in s=m^2, location m0^2, scale m0*Gamma0), truncated
// to +-3*Gamma in mass, via the exact inverse CDF in s. Consumes EXACTLY ONE counter-RNG draw
// for ANY meson (even pseudoscalars, Gamma=0) so the GPU/CPU RNG streams stay phase-aligned.
#pragma once
#include <cmath>
#include <cstdint>
#include "../common/rng.cuh"

__host__ __device__ inline double mesonWidth(int pdg){
  switch(abs(pdg)){
    case 213: case 113: return 0.1491;   // rho+-, rho0
    case 223: return 0.00849;            // omega
    case 333: return 0.00425;            // phi
    case 323: return 0.0473;             // K*+-
    case 313: return 0.0487;             // K*0
  } return 0.0;                           // pseudoscalars: zero width -> pole mass
}
// Sample the meson mass: BW for vectors (one draw), pole for pseudoscalars (draw consumed anyway).
__host__ __device__ inline double sampleBWmass(int pdg,double m0,uint64_t& ctr){
  double r=u01(splitmix64(ctr++));                 // ALWAYS one draw -> stream stays aligned
  double G=mesonWidth(pdg);
  if(G<=0.0) return m0;
  const double N=3.0;
  double lo=m0-N*G; if(lo<0.0) lo=0.0;
  double sLo=lo*lo, sHi=(m0+N*G)*(m0+N*G), s0=m0*m0, gw=m0*G;
  double aLo=atan((sLo-s0)/gw), aHi=atan((sHi-s0)/gw);
  return sqrt(s0 + gw*tan(aLo + r*(aHi-aLo)));
}
