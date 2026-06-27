// cuPythia — GPU recursive hadron decays (behind -DDECAYS). Decays the primary unstable hadrons
// the Lund hadronizer makes (rho/K*/omega/phi/eta/eta'/K0 -> pi/K/gamma/K_S/K_L) into the ALEPH
// particle-level stable set, one event per thread, recursion-free (LIFO stack), counter-RNG so
// host==device and events stay O(1)-reproducible, 4-momentum conserved per decay by construction.
//
// This closes the #1 gap vs real ALEPH LEP1 data (charged multiplicity 11.3 -> ~20): without decays
// rho->pipi, K*->Kpi, omega->3pi etc. never produce the extra (charged) tracks a detector sees.
//
// SCOPE (honest, v1): light-unstable-meson 2-/3-body FLAT phase space; BRs from Pythia 8.317
// ParticleData.xml, truncated to kept channels + renormalized to 1.0. pi0 and K_S kept STABLE
// (Rivet particle-level / ALEPH unfolded convention). No D/B (declared stable), no baryons, no
// Dalitz/angular ME shapes, no rare leptonic tails -> a documented few-% residual.
#pragma once
#include <cmath>
#include <cstdint>
#include "../common/rng.cuh"
#include "bw_inc.cuh"       // sampleBWmass (always 1 draw), mesonWidth
#include "shower_inc.cuh"   // boostBy, rotZto

#define MAXFINAL      256
#define MAXSTACK      96
#define MAXDECAYITERS 512
#define KDAL          8

// Pole masses for decay products/parents — MUST match mesonMassMR (hadronize_mr.cu) for the shared
// ids so the kern's on-shell check passes; plus gamma/K_L/K_S.
__host__ __device__ inline double dPoleMass(int ap){
  switch(ap){
    case 22:  return 0.0;
    case 111: return 0.13498; case 211: return 0.13957;
    case 321: return 0.49368; case 311: return 0.49761; case 130: return 0.49761; case 310: return 0.49761;
    case 221: return 0.54786; case 331: return 0.95778;
    case 113: return 0.77526; case 213: return 0.77526; case 223: return 0.78266; case 333: return 1.01946;
    case 313: return 0.89555; case 323: return 0.89167;
  } return 0.0;
}
// Charge conjugate (flip non-self-conjugate ids; keep pi0/eta/eta'/rho0/omega/phi/gamma/K_L/K_S).
__host__ __device__ inline int ccDec(int p){
  switch(p){ case 211:return -211; case -211:return 211; case 321:return -321; case -321:return 321;
    case 213:return -213; case -213:return 213; case 323:return -323; case -323:return 323;
    case 311:return -311; case -311:return 311; case 313:return -313; case -313:return 313; }
  return p;
}

// First-pass decay table, encoded as switch/local-constexpr accessors (file-scope struct/arrays are
// NOT visible to device code under nvcc; function-local constexpr scalar tables are usable on BOTH
// host and device). 9 parent rows, 32 channels. Per-channel BR renormalized to 1.0 per parent;
// products SIGNED PDG; antiparticle parents (pid<0) conjugate products at runtime via ccDec.
// Rows: 0 eta(221) 1 eta'(331) 2 rho0(113) 3 rho+(213) 4 omega(223) 5 phi(333) 6 K*0(313) 7 K*+(323) 8 K0(311).
__host__ __device__ inline int dRow(int pdg){
  constexpr int PP[9]={221,331,113,213,223,333,313,323,311};
  int a=abs(pdg); for(int i=0;i<9;++i) if(PP[i]==a) return i; return -1;   // -1 = STABLE
}
__host__ __device__ inline void dParentInfo(int row,int& first,int& nch){
  constexpr int PF[9]={0,4,9,12,14,17,24,27,30}, PN[9]={4,5,3,2,3,7,3,3,2};
  first=PF[row]; nch=PN[row];
}
__host__ __device__ inline double dChanBR(int ci){
  constexpr double B[32]={0.3931,0.3257,0.2274,0.0460, 0.4366,0.2947,0.2173,0.0277,0.0219,
    0.9988,0.0006,0.0006, 0.99955,0.00045, 0.8995,0.0835,0.0154,
    0.4893,0.3422,0.0421,0.0421,0.0421,0.0270,0.0131, 0.6649,0.3327,0.0024, 0.6660,0.3330,0.0010, 0.5,0.5};
  return B[ci];
}
__host__ __device__ inline int dChanN(int ci){
  constexpr int NP[32]={2,3,3,3, 3,2,3,2,2, 2,2,2, 2,2, 3,2,2, 2,2,2,2,2,3,2, 2,2,2, 2,2,2, 1,1};
  return NP[ci];
}
__host__ __device__ inline int dChanProd(int ci,int j){
  constexpr int P0[32]={22,111,211,211, 211,113,111,223,22, 211,111,221, 211,211, 211,111,211,
    321,130,-213,113,213,211,221, 321,311,311, 311,321,321, 130,310};
  constexpr int P1[32]={22,111,-211,-211, -211,22,111,22,22, -211,22,22, 111,22, -211,22,-211,
    -321,310,211,111,-211,-211,22, -211,111,22, 211,111,22, 0,0};
  constexpr int P2[32]={0,111,111,22, 221,0,221,0,0, 0,0,0, 0,0, 111,0,0, 0,0,0,0,0,111,0, 0,0,0, 0,0,0, 0,0};
  return (j==0)?P0[ci]:((j==1)?P1[ci]:P2[ci]);
}
// Mass of a product (BW for vectors via bw_inc, pole otherwise) — ALWAYS exactly one RNG draw.
__host__ __device__ inline double decayMass(int pid,uint64_t& ctr){ return sampleBWmass(pid, dPoleMass(abs(pid)), ctr); }

// Isotropic 2-body in the parent rest frame (sum p=0, sum E=M). Exactly 2 draws (cos,phi).
__host__ __device__ inline bool twoBody(double M,double m1,double m2,uint64_t& ctr,double* q1,double* q2){
  double cth=2.0*u01(splitmix64(ctr++))-1.0, ph=2.0*M_PI*u01(splitmix64(ctr++));   // ALWAYS 2 draws
  if(M<m1+m2) return false;
  double lam=(M*M-m1*m1-m2*m2)*(M*M-m1*m1-m2*m2)-4.0*m1*m1*m2*m2;
  double pst=0.5*sqrt(fmax(0.0,lam))/M, E1=(M*M+m1*m1-m2*m2)/(2.0*M), E2=M-E1;
  double sth=sqrt(fmax(0.0,1.0-cth*cth));
  q1[0]=pst*sth*cos(ph); q1[1]=pst*sth*sin(ph); q1[2]=pst*cth; q1[3]=E1;
  q2[0]=-q1[0]; q2[1]=-q1[1]; q2[2]=-q1[2]; q2[3]=E2;
  return true;
}
// Dalitz physical region test (s12=(p1+p2)^2, s23=(p2+p3)^2).
__host__ __device__ inline bool inDalitz(double M,double m1,double m2,double m3,double s12,double s23){
  if(s12<(m1+m2)*(m1+m2) || s12>(M-m3)*(M-m3)) return false;
  double rs=sqrt(s12), E2=(s12-m1*m1+m2*m2)/(2.0*rs), E3=(M*M-s12-m3*m3)/(2.0*rs);
  double p2=E2*E2-m2*m2, p3=E3*E3-m3*m3; if(p2<0||p3<0) return false; p2=sqrt(p2); p3=sqrt(p3);
  double smin=(E2+E3)*(E2+E3)-(p2+p3)*(p2+p3), smax=(E2+E3)*(E2+E3)-(p2-p3)*(p2-p3);
  return s23>=smin && s23<=smax;
}
// Flat 3-body in the parent rest frame, random isotropic orientation. FIXED 2*KDAL+2 draws (all KDAL
// trials executed -> phase invariant; keep first valid Dalitz point; +2 for random orientation).
__host__ __device__ inline bool threeBody(double M,double m1,double m2,double m3,uint64_t& ctr,
                                          double* q1,double* q2,double* q3){
  double s12lo=(m1+m2)*(m1+m2), s12hi=(M-m3)*(M-m3), s23lo=(m2+m3)*(m2+m3), s23hi=(M-m1)*(M-m1);
  bool found=false; double S12=0,S23=0;
  for(int t=0;t<KDAL;++t){ double a=u01(splitmix64(ctr++)), b=u01(splitmix64(ctr++));
    double s12=s12lo+a*(s12hi-s12lo), s23=s23lo+b*(s23hi-s23lo);
    if(!found && M>m1+m2+m3 && inDalitz(M,m1,m2,m3,s12,s23)){ S12=s12; S23=s23; found=true; } }
  double cth=2.0*u01(splitmix64(ctr++))-1.0, ph=2.0*M_PI*u01(splitmix64(ctr++));   // ALWAYS 2 (orientation)
  if(!found) return false;
  double M2=M*M, E1=(M2+m1*m1-S23)/(2.0*M), E3=(M2+m3*m3-S12)/(2.0*M), E2=M-E1-E3;
  double p1=sqrt(fmax(0.0,E1*E1-m1*m1)), p2=sqrt(fmax(0.0,E2*E2-m2*m2)), p3=sqrt(fmax(0.0,E3*E3-m3*m3));
  double c12=(p1>0&&p2>0)?(p3*p3-p1*p1-p2*p2)/(2.0*p1*p2):0.0; if(c12>1)c12=1; if(c12<-1)c12=-1;
  double s12a=sqrt(fmax(0.0,1.0-c12*c12));
  q1[0]=0;        q1[1]=0; q1[2]=p1;       q1[3]=E1;
  q2[0]=p2*s12a;  q2[1]=0; q2[2]=p2*c12;   q2[3]=E2;
  q3[0]=-(q1[0]+q2[0]); q3[1]=0; q3[2]=-(q1[2]+q2[2]); q3[3]=E3;
  double sth=sqrt(fmax(0.0,1.0-cth*cth)), nx=sth*cos(ph), ny=sth*sin(ph), nz=cth;
  rotZto(nx,ny,nz,q1); rotZto(nx,ny,nz,q2); rotZto(nx,ny,nz,q3);   // random isotropic orientation
  return true;
}

// Decay one event's primary hadrons (H/hid, nH) into the final stable list (F/fid). Returns nF, or
// -1 on overflow / kinematic failure (caller drops the event). LIFO stack, no recursion, no alloc.
__host__ __device__ inline int decayEvent(const double* H,const int* hid,int nH,uint64_t ctr,
                                          double* F,int* fid){
  double S[MAXSTACK*4]; int Sid[MAXSTACK]; int sp=0;
  for(int i=nH-1;i>=0;--i){ if(sp>=MAXSTACK) return -1; for(int k=0;k<4;++k)S[4*sp+k]=H[4*i+k]; Sid[sp]=hid[i]; sp++; }
  int nF=0, guard=0;
  while(sp>0 && guard<MAXDECAYITERS){ guard++;
    sp--; double P[4]={S[4*sp],S[4*sp+1],S[4*sp+2],S[4*sp+3]}; int pid=Sid[sp];
    int row=dRow(pid);
    if(row<0){ if(nF>=MAXFINAL) return -1; for(int k=0;k<4;++k)F[4*nF+k]=P[k]; fid[nF]=pid; nF++; continue; }
    double r=u01(splitmix64(ctr++));                                  // channel pick: 1 draw
    int first,nc; dParentInfo(row,first,nc); int ci=first+nc-1; double cum=0;
    for(int c=0;c<nc;++c){ cum+=dChanBR(first+c); if(r<cum){ ci=first+c; break; } }
    int np=dChanN(ci), prod[3]; double m[3];
    for(int j=0;j<np;++j){ int pp=dChanProd(ci,j); if(pid<0) pp=ccDec(pp); prod[j]=pp; m[j]=decayMass(pp,ctr); } // np draws
    double M=sqrt(fmax(0.0,P[3]*P[3]-P[0]*P[0]-P[1]*P[1]-P[2]*P[2]));
    double q[3][4]; bool ok=true;
    if(np==1){ for(int k=0;k<4;++k)q[0][k]=P[k]; }                    // relabel (K0->K_S/K_L): copy 4-vec
    else if(np==2){ ok=twoBody(M,m[0],m[1],ctr,q[0],q[1]); }
    else          { ok=threeBody(M,m[0],m[1],m[2],ctr,q[0],q[1],q[2]); }
    if(!ok) return -1;
    if(np>1){ double bx=P[0]/P[3],by=P[1]/P[3],bz=P[2]/P[3],ga=P[3]/M;   // rest frame -> lab
      for(int j=0;j<np;++j){ double o[4]; boostBy(q[j],bx,by,bz,ga,o); for(int k=0;k<4;++k)q[j][k]=o[k]; } }
    for(int j=0;j<np;++j){ if(sp>=MAXSTACK) return -1; for(int k=0;k<4;++k)S[4*sp+k]=q[j][k]; Sid[sp]=prod[j]; sp++; }
  }
  return (guard>=MAXDECAYITERS)? -1 : nF;
}
