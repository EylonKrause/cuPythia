// cuPythia — GPU recursive hadron decays (behind -DDECAYS). Decays the primary unstable hadrons
// the Lund hadronizer makes (rho/K*/omega/phi/eta/eta'/K0 -> pi/K/gamma/K_S/K_L) into the ALEPH
// particle-level stable set, one event per thread, recursion-free (LIFO stack), counter-RNG so
// host==device and events stay O(1)-reproducible, 4-momentum conserved per decay by construction.
//
// This closes the #1 gap vs real ALEPH LEP1 data (charged multiplicity 11.3 -> ~20): without decays
// rho->pipi, K*->Kpi, omega->3pi etc. never produce the extra (charged) tracks a detector sees.
//
// SCOPE (honest): light-unstable-meson 2-/3-body FLAT phase space; BRs from Pythia 8.317
// ParticleData.xml, truncated to kept channels + renormalized to 1.0. pi0 and K_S kept STABLE
// (Rivet particle-level / ALEPH unfolded convention). Optional extensions, each opt-in so the
// -DDECAYS-only build stays byte-identical:
//   -DBARYONS    baryon-resonance decay rows (Delta/Sigma*/Xi*/Sigma0).
//   -DHFDECAY    D/B heavy-flavour decays (separate table + fixed-budget fourBody sampler): D*/Ds*
//                radiative+pi, D/Ds truncated weak Cabibbo set, effective B->D(*)+npi (capped 4-body,
//                <n_ch>_B ~4.1 vs PDG 4.97 -- a documented undershoot).
//   -DDALITZ_ME  real Dalitz densities for omega/phi (P-wave |p+ x p-|^2) and eta (linear slope) via
//                accept-reject, instead of flat 3-body.
// Still out of scope: rare leptonic tails, charm/bottom BARYONS, 5-body+ B channels -> a few-% residual.
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
    // open charm/bottom (parents + intermediate products) -- only reached under -DHFDECAY
    case 411: return 1.86966; case 421: return 1.86484; case 431: return 1.96835;
    case 413: return 2.01026; case 423: return 2.00685; case 433: return 2.11220;
    case 511: return 5.27966; case 521: return 5.27934; case 531: return 5.36688;
    case 513: return 5.32471; case 523: return 5.32471; case 533: return 5.41550;
    // stable baryons (decay products) -- only reached under -DBARYONS
    case 2212: return 0.93827; case 2112: return 0.93957; case 3122: return 1.11568;
    case 3112: return 1.19745; case 3212: return 1.19264; case 3222: return 1.18937;
    case 3312: return 1.32171; case 3322: return 1.31486; case 3334: return 1.67245;
  } return 0.0;
}
// Charge conjugate (flip non-self-conjugate ids; keep pi0/eta/eta'/rho0/omega/phi/gamma/K_L/K_S).
__host__ __device__ inline int ccDec(int p){
  switch(p){ case 211:return -211; case -211:return 211; case 321:return -321; case -321:return 321;
    case 213:return -213; case -213:return 213; case 323:return -323; case -323:return 323;
    case 311:return -311; case -311:return 311; case 313:return -313; case -313:return 313;
    // open charm/bottom (B/D* cascade products + B* parents) -- only reached under -DHFDECAY
    case 411:return -411; case -411:return 411; case 421:return -421; case -421:return 421;
    case 431:return -431; case -431:return 431; case 413:return -413; case -413:return 413;
    case 423:return -423; case -423:return 423; case 433:return -433; case -433:return 433;
    case 511:return -511; case -511:return 511; case 521:return -521; case -521:return 521;
    case 531:return -531; case -531:return 531;
    // baryon products (antibaryon parents conjugate their products) -- only under -DBARYONS
    case 2212:return -2212; case -2212:return 2212; case 2112:return -2112; case -2112:return 2112;
    case 3122:return -3122; case -3122:return 3122; case 3222:return -3222; case -3222:return 3222;
    case 3212:return -3212; case -3212:return 3212; case 3112:return -3112; case -3112:return 3112;
    case 3312:return -3312; case -3312:return 3312; case 3322:return -3322; case -3322:return 3322; }
  return p;
}

// First-pass decay table, encoded as switch/local-constexpr accessors (file-scope struct/arrays are
// NOT visible to device code under nvcc; function-local constexpr scalar tables are usable on BOTH
// host and device). 9 parent rows, 32 channels. Per-channel BR renormalized to 1.0 per parent;
// products SIGNED PDG; antiparticle parents (pid<0) conjugate products at runtime via ccDec.
// Rows: 0 eta(221) 1 eta'(331) 2 rho0(113) 3 rho+(213) 4 omega(223) 5 phi(333) 6 K*0(313) 7 K*+(323) 8 K0(311).
// Under -DBARYONS the table gains 10 baryon-resonance parents (rows 9-18) / 20 channels (32-51):
// Delta++,Delta+,Delta0,Delta- -> N pi; Sigma*+/0/- -> Lambda/Sigma pi; Xi*0/- -> Xi pi; Sigma0 -> Lambda gamma.
// (The weak baryons p,n,Lambda,Sigma+-,Xi-,Xi0,Omega- are kept STABLE at particle level, no rows.)
__host__ __device__ inline int dRow(int pdg){
  int a=abs(pdg);
#ifdef BARYONS
  constexpr int PP[19]={221,331,113,213,223,333,313,323,311, 2224,2214,2114,1114,3224,3214,3114,3324,3314,3212};
  for(int i=0;i<19;++i) if(PP[i]==a) return i; return -1;
#else
  constexpr int PP[9]={221,331,113,213,223,333,313,323,311};
  for(int i=0;i<9;++i) if(PP[i]==a) return i; return -1;   // -1 = STABLE
#endif
}
__host__ __device__ inline void dParentInfo(int row,int& first,int& nch){
#ifdef BARYONS
  constexpr int PF[19]={0,4,9,12,14,17,24,27,30, 32,33,35,37,38,41,44,47,49,51};
  constexpr int PN[19]={4,5,3,2,3,7,3,3,2, 1,2,2,1,3,3,3,2,2,1};
#else
  constexpr int PF[9]={0,4,9,12,14,17,24,27,30}, PN[9]={4,5,3,2,3,7,3,3,2};
#endif
  first=PF[row]; nch=PN[row];
}
__host__ __device__ inline double dChanBR(int ci){
  constexpr double B[52]={0.3931,0.3257,0.2274,0.0460, 0.4366,0.2947,0.2173,0.0277,0.0219,
    0.9988,0.0006,0.0006, 0.99955,0.00045, 0.8995,0.0835,0.0154,
    0.4893,0.3422,0.0421,0.0421,0.0421,0.0270,0.0131, 0.6649,0.3327,0.0024, 0.6660,0.3330,0.0010, 0.5,0.5,
    1.0, 0.667,0.333, 0.667,0.333, 1.0, 0.88,0.06,0.06, 0.88,0.06,0.06, 0.88,0.06,0.06, 0.67,0.33, 0.67,0.33, 1.0};
  return B[ci];
}
__host__ __device__ inline int dChanN(int ci){
  constexpr int NP[52]={2,3,3,3, 3,2,3,2,2, 2,2,2, 2,2, 3,2,2, 2,2,2,2,2,3,2, 2,2,2, 2,2,2, 1,1,
    2, 2,2, 2,2, 2, 2,2,2, 2,2,2, 2,2,2, 2,2, 2,2, 2};
  return NP[ci];
}
__host__ __device__ inline int dChanProd(int ci,int j){
  constexpr int P0[52]={22,111,211,211, 211,113,111,223,22, 211,111,221, 211,211, 211,111,211,
    321,130,-213,113,213,211,221, 321,311,311, 311,321,321, 130,310,
    2212, 2212,2112, 2112,2212, 2112, 3122,3222,3212, 3122,3222,3112, 3122,3212,3112, 3312,3322, 3322,3312, 3122};
  constexpr int P1[52]={22,111,-211,-211, -211,22,111,22,22, -211,22,22, 111,22, -211,22,-211,
    -321,310,211,111,-211,-211,22, -211,111,22, 211,111,22, 0,0,
    211, 111,211, 111,-211, -211, 211,111,211, 111,-211,211, -211,-211,111, 211,111, -211,111, 22};
  constexpr int P2[52]={0,111,111,22, 221,0,221,0,0, 0,0,0, 0,0, 111,0,0, 0,0,0,0,0,111,0, 0,0,0, 0,0,0, 0,0,
    0, 0,0, 0,0, 0, 0,0,0, 0,0,0, 0,0,0, 0,0, 0,0, 0};
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
#ifdef DALITZ_ME
#define KDAL_ME 32     // more Dalitz trials when ME-weighting (accept-reject) so flat-fallback is rare
// Normalized Dalitz matrix-element weight in [0,1] for a candidate (s12,s23). shape 0 = flat (1.0);
// shape 1 = omega/phi -> pi+pi-pi0 P-wave |p1 x p2|^2 (q1,q2 charged); shape 2 = eta -> pi+pi-pi0
// linear+quadratic density 1 + a y + b y^2 + d x^2 (PDG charged-channel a,b,d). Envelopes are rigorous
// upper bounds (clamped), so accept-reject is unbiased; all-reject falls back to flat (now <0.01%).
__host__ __device__ inline double dalitzW(int shape,double M,double m1,double m2,double m3,double s12,double s23){
  if(shape==0) return 1.0;
  double M2=M*M, E1=(M2+m1*m1-s23)/(2.0*M), E3=(M2+m3*m3-s12)/(2.0*M), E2=M-E1-E3;
  double p1=sqrt(fmax(0.0,E1*E1-m1*m1)), p2=sqrt(fmax(0.0,E2*E2-m2*m2)), p3=sqrt(fmax(0.0,E3*E3-m3*m3));
  if(shape==1){
    double dot=0.5*(p3*p3-p1*p1-p2*p2);                              // p1.p2 (rest frame: p1+p2+p3=0)
    double cr=p1*p1*p2*p2-dot*dot; if(cr<0.0) cr=0.0;                // |p1 x p2|^2
    double u23=m2+m3, u13=m1+m3;                                     // envelope base: (p1max p2max)^2
    double l1=(M2-(m1+u23)*(m1+u23))*(M2-(m1-u23)*(m1-u23));
    double l2=(M2-(m2+u13)*(m2+u13))*(M2-(m2-u13)*(m2-u13));
    double p1m2=0.25*fmax(0.0,l1)/M2, p2m2=0.25*fmax(0.0,l2)/M2;
    // p1 and p2 cannot both be maximal at once, so (p1max p2max)^2 over-counts ~6x. The numerically-
    // verified GLOBAL max of cr/(p1max^2 p2max^2) over the Dalitz region is 0.157 (omega) / 0.154 (phi)
    // / 0.148 (massless limit); the 0.20 factor is a rigorous tightened envelope (>0.16) that keeps the
    // clamp from firing (unbiased) yet raises acceptance ~5x -> flat-fallback drops from ~20% to <0.01%.
    double wmax=0.20*p1m2*p2m2;
    if(wmax<=0.0) return 0.0; double w=cr/wmax; return w>1.0?1.0:w;
  }
  double Q=M-m1-m2-m3; if(Q<=0.0) return 1.0;
  double T1=E1-m1, T2=E2-m2, T3=E3-m3;                               // kinetic energies (q3 = pi0)
  double y=3.0*T3/Q-1.0, x=1.7320508075688772*(T1-T2)/Q;            // sqrt(3)
  const double a=-1.095, b=0.145, d=0.081;                          // PDG eta->pi+pi-pi0 Dalitz params
  double w=1.0+a*y+b*y*y+d*x*x; if(w<0.0) w=0.0;
  double wmax=1.0+1.095+0.145+0.081+0.10;                           // bound at y=-1,|x|=1 (+0.10 safety)
  double r=w/wmax; return r>1.0?1.0:r;
}
// Pick the Dalitz shape for a 3-body channel: omega(223)/phi(333)->pi+pi-pi0 = 1; eta(221)->pi+pi-pi0 = 2.
__host__ __device__ inline int dalitzShapeFor(int pid,int p0,int p1,int p2){
  int a=abs(pid); int n211=(abs(p0)==211)+(abs(p1)==211)+(abs(p2)==211);
  int n111=(abs(p0)==111)+(abs(p1)==111)+(abs(p2)==111);
  if(n211==2&&n111==1){ if(a==223||a==333) return 1; if(a==221) return 2; }
  return 0;
}
#endif
// 3-body in the parent rest frame, random isotropic orientation. FIXED draws (all trials executed ->
// phase invariant; keep first valid Dalitz point; +2 for orientation). Default: flat phase space
// (2*KDAL+2 draws). Under -DDALITZ_ME: accept-reject the ME shape (3*KDAL_ME+2 draws); shape==0 -> flat.
__host__ __device__ inline bool threeBody(double M,double m1,double m2,double m3,uint64_t& ctr,
                                          double* q1,double* q2,double* q3,int shape=0){
  double s12lo=(m1+m2)*(m1+m2), s12hi=(M-m3)*(M-m3), s23lo=(m2+m3)*(m2+m3), s23hi=(M-m1)*(M-m1);
  bool found=false; double S12=0,S23=0;
#ifdef DALITZ_ME
  bool foundME=false; double Me12=0,Me23=0; const int NTRY=KDAL_ME;
#else
  const int NTRY=KDAL; (void)shape;
#endif
  for(int t=0;t<NTRY;++t){ double a=u01(splitmix64(ctr++)), b=u01(splitmix64(ctr++));
    double s12=s12lo+a*(s12hi-s12lo), s23=s23lo+b*(s23hi-s23lo);
    bool valid=(M>m1+m2+m3)&&inDalitz(M,m1,m2,m3,s12,s23);
    if(!found && valid){ S12=s12; S23=s23; found=true; }
#ifdef DALITZ_ME
    double w=valid?dalitzW(shape,M,m1,m2,m3,s12,s23):0.0;
    double ra=u01(splitmix64(ctr++));                               // accept draw (FIXED 3rd draw/trial)
    if(!foundME && valid && ra<w){ Me12=s12; Me23=s23; foundME=true; }
#endif
  }
#ifdef DALITZ_ME
  if(foundME){ S12=Me12; S23=Me23; }                                // prefer ME-accepted; else flat fallback
#endif
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
#ifdef HFDECAY
// Sequential (iterated 2-body) 4-body phase space: M -> (123)+4, (123) -> (12)+3, (12) -> 1+2.
// Intermediate masses by nested-uniform sampling (always in the physical range), each stage an
// isotropic 2-body boosted into the lab. FIXED 8 draws (2 masses + 3 x twoBody). NOT flat 4-body
// phase space (it is a valid, 4-momentum-exact kinematic model adequate for track counting; the
// flat-4-body Dalitz shape is out of scope) -- used only for the effective D/B multi-pion channels.
__host__ __device__ inline bool fourBody(double M,double m1,double m2,double m3,double m4,uint64_t& ctr,
                                         double* q1,double* q2,double* q3,double* q4){
  double r1=u01(splitmix64(ctr++)), r2=u01(splitmix64(ctr++));            // 2 mass draws
  bool ok=(M>m1+m2+m3+m4);
  double m12 =ok?((m1+m2)+r1*((M-m3-m4)-(m1+m2))):1.0;
  double m123=ok?((m12+m3)+r2*((M-m4)-(m12+m3))):1.0;
  double a[4],b[4],c[4],d[4],e[4],f[4];
  bool s1=twoBody(ok?M   :1.0, ok?m123:0.0, ok?m4:0.0, ctr,a,b);          // M -> (123)+4   (2 draws)
  bool s2=twoBody(ok?m123:1.0, ok?m12 :0.0, ok?m3:0.0, ctr,c,d);          // (123)->(12)+3  (2 draws)
  bool s3=twoBody(ok?m12 :1.0, ok?m1  :0.0, ok?m2:0.0, ctr,e,f);          // (12)->1+2      (2 draws)
  if(!(ok&&s1&&s2&&s3)) return false;
  for(int k=0;k<4;++k) q4[k]=b[k];
  double bx=a[0]/a[3],by=a[1]/a[3],bz=a[2]/a[3],ga=a[3]/m123, cl[4],dl[4];
  boostBy(c,bx,by,bz,ga,cl); boostBy(d,bx,by,bz,ga,dl); for(int k=0;k<4;++k) q3[k]=dl[k];
  double bx2=cl[0]/cl[3],by2=cl[1]/cl[3],bz2=cl[2]/cl[3],ga2=cl[3]/m12;
  boostBy(e,bx2,by2,bz2,ga2,q1); boostBy(f,bx2,by2,bz2,ga2,q2);
  return true;
}
// ---- D/B decay table (separate accessors so -DDECAYS-only stays byte-identical). 12 parents, 30
// channels: D*/Ds* radiative+pi (2-body); D0/D+/Ds weak (truncated+renormalized 2-/3-/4-body Cabibbo-
// favoured set, BRs from PDG); B0/B+/Bs EFFECTIVE B->D(*)+npi tuned toward <n_ch>_B (capped at 4-body
// -> ~4.1 vs PDG 4.97, a documented ~15% undershoot). Products SIGNED for the POSITIVE parent; pid<0
// conjugated by ccDec. Rows: 0 D*+(413) 1 D*0(423) 2 D0(421) 3 D+(411) 4 Ds+(431) 5 Ds*+(433)
//        6 B0(511) 7 B+(521) 8 Bs0(531) 9 B*0(513) 10 B*+(523) 11 Bs*0(533).
__host__ __device__ inline int hfRow(int pdg){ int a=abs(pdg);
  constexpr int PP[12]={413,423,421,411,431,433, 511,521,531,513,523,533};
  for(int i=0;i<12;++i) if(PP[i]==a) return i; return -1; }
__host__ __device__ inline void hfParentInfo(int row,int& first,int& nch){
  constexpr int PF[12]={0,3,5,8,11,14, 16,20,24,27,28,29};
  constexpr int PN[12]={3,2,3,3,3,2, 4,4,3,1,1,1};
  first=PF[row]; nch=PN[row]; }
__host__ __device__ inline double hfChanBR(int ci){
  constexpr double B[30]={0.677,0.307,0.016, 0.647,0.353, 0.1488,0.5424,0.3088, 0.5054,0.3333,0.1613,
    0.5745,0.1170,0.3085, 0.9416,0.0584,
    0.30,0.25,0.20,0.25, 0.30,0.25,0.20,0.25, 0.40,0.30,0.30, 1.0, 1.0, 1.0};
  return B[ci]; }
__host__ __device__ inline int hfChanN(int ci){
  constexpr int NP[30]={2,2,2, 2,2, 2,3,4, 3,4,2, 3,3,2, 2,2,
    4,2,2,3, 4,2,2,3, 4,2,3, 2, 2, 2};
  return NP[ci]; }
__host__ __device__ inline int hfChanProd(int ci,int j){
  constexpr int P0[30]={421,411,411, 421,421, -321,-321,-321, -321,-321,-311, 321,211,321, 431,431,
    -411,-413,-411,-411, -421,-423,-421,-421, -431,-431,-431, 511, 521, 531};
  constexpr int P1[30]={211,111,22, 111,22, 211,211,211, 211,211,211, -321,211,-311, 22,111,
    211,211,211,211, 211,211,211,211, 211,211,211, 22, 22, 22};
  constexpr int P2[30]={0,0,0, 0,0, 0,111,211, 211,211,0, 211,-211,0, 0,0,
    211,0,0,111, 211,0,0,111, 211,0,111, 0, 0, 0};
  constexpr int P3[30]={0,0,0, 0,0, 0,0,-211, 0,111,0, 0,0,0, 0,0,
    -211,0,0,0, -211,0,0,0, -211,0,0, 0, 0, 0};
  return (j==0)?P0[ci]:((j==1)?P1[ci]:((j==2)?P2[ci]:P3[ci])); }
#endif

// Decay one event's primary hadrons (H/hid, nH) into the final stable list (F/fid). Returns nF, or
// -1 on overflow / kinematic failure (caller drops the event). LIFO stack, no recursion, no alloc.
__host__ __device__ __noinline__ int decayEvent(const double* H,const int* hid,int nH,uint64_t ctr,
                                          double* F,int* fid){
  double S[MAXSTACK*4]; int Sid[MAXSTACK]; int sp=0;
  for(int i=nH-1;i>=0;--i){ if(sp>=MAXSTACK) return -1; for(int k=0;k<4;++k)S[4*sp+k]=H[4*i+k]; Sid[sp]=hid[i]; sp++; }
  int nF=0, guard=0;
  while(sp>0 && guard<MAXDECAYITERS){ guard++;
    sp--; double P[4]={S[4*sp],S[4*sp+1],S[4*sp+2],S[4*sp+3]}; int pid=Sid[sp];
    int np=0, prod[4]; double m[4]; bool isParent=false;
#ifdef HFDECAY
    int hrow=hfRow(pid);                                              // D/B parents (separate table) first
    if(hrow>=0){ isParent=true;
      double r=u01(splitmix64(ctr++));                               // channel pick: 1 draw
      int first,nc; hfParentInfo(hrow,first,nc); int ci=first+nc-1; double cum=0;
      for(int c=0;c<nc;++c){ cum+=hfChanBR(first+c); if(r<cum){ ci=first+c; break; } }
      np=hfChanN(ci);
      for(int j=0;j<np;++j){ int pp=hfChanProd(ci,j); if(pid<0) pp=ccDec(pp); prod[j]=pp; m[j]=decayMass(pp,ctr); }
    } else
#endif
    { int row=dRow(pid);
      if(row<0){ if(nF>=MAXFINAL) return -1; for(int k=0;k<4;++k)F[4*nF+k]=P[k]; fid[nF]=pid; nF++; continue; }
      isParent=true;
      double r=u01(splitmix64(ctr++));                                // channel pick: 1 draw
      int first,nc; dParentInfo(row,first,nc); int ci=first+nc-1; double cum=0;
      for(int c=0;c<nc;++c){ cum+=dChanBR(first+c); if(r<cum){ ci=first+c; break; } }
      np=dChanN(ci);
      for(int j=0;j<np;++j){ int pp=dChanProd(ci,j); if(pid<0) pp=ccDec(pp); prod[j]=pp; m[j]=decayMass(pp,ctr); } // np draws
    }
    (void)isParent;
    double M=sqrt(fmax(0.0,P[3]*P[3]-P[0]*P[0]-P[1]*P[1]-P[2]*P[2]));
    double q[4][4]; bool ok=true;
    if(np==1){ for(int k=0;k<4;++k)q[0][k]=P[k]; }                    // relabel (K0->K_S/K_L): copy 4-vec
    else if(np==2){ ok=twoBody(M,m[0],m[1],ctr,q[0],q[1]); }
#ifdef DALITZ_ME
    else if(np==3){ ok=threeBody(M,m[0],m[1],m[2],ctr,q[0],q[1],q[2],dalitzShapeFor(pid,prod[0],prod[1],prod[2])); }
#else
    else if(np==3){ ok=threeBody(M,m[0],m[1],m[2],ctr,q[0],q[1],q[2]); }
#endif
#ifdef HFDECAY
    else          { ok=fourBody(M,m[0],m[1],m[2],m[3],ctr,q[0],q[1],q[2],q[3]); }
#else
    else          { ok=false; }
#endif
    if(!ok) return -1;
    if(np>1){ double bx=P[0]/P[3],by=P[1]/P[3],bz=P[2]/P[3],ga=P[3]/M;   // rest frame -> lab
      for(int j=0;j<np;++j){ double o[4]; boostBy(q[j],bx,by,bz,ga,o); for(int k=0;k<4;++k)q[j][k]=o[k]; } }
    for(int j=0;j<np;++j){ if(sp>=MAXSTACK) return -1; for(int k=0;k<4;++k)S[4*sp+k]=q[j][k]; Sid[sp]=prod[j]; sp++; }
  }
  return (guard>=MAXDECAYITERS)? -1 : nF;
}
