// Shared FSR dipole-shower core (extracted verbatim from shower_fsr.cu so the shower and the
// shower->hadronization bridge use the SAME validated code). One event per thread, counter-RNG,
// large-N_c colour chain (partons kept in colour order q ... gluons ... qbar). See shower_fsr.cu
// for the full validation (thrust vs Pythia LL 4.0%, 100% GPU==CPU control flow).
#pragma once
#include <cmath>
#include <cstdint>
#include "../common/rng.cuh"

static const double MZ     = 91.1876;       // Z mass / CM energy
static const double EBEAM  = 0.5*MZ;        // each initial quark energy
static const double PT2MIN = 0.5*0.5;       // (TimeShower:pTmin = 0.5 GeV)^2 cutoff
static const double ASMZ   = 0.1365;        // alpha_s(M_Z)
#define MAXP 64
#ifdef G2QQ_WEIGHT4
#define G2Q4_OVER 4.0   // veto oversample for the option-4 g->qqbar weight (zCosThe makes it >1 in places)
#endif

// Running alpha_s, flavour-threshold matched (n_f=5,4,3 across m_b,m_c). 1-loop by default;
// -DAS_2LOOP selects the 2-loop form matching Pythia AlphaStrong order=2 (bit-identical, see
// as2loop_validate.cc). NOTE (honest): a leading-log shower conventionally uses 1-loop(+CMW);
// 2-loop is offered for completeness, not because it is the right order for THIS shower.
__host__ __device__ inline double alphaS(double mu2){
  const double mc2=1.5*1.5, mb2=4.8*4.8;
#ifdef AS_2LOOP
  // 2-loop alpha_s (RPP 2006 eq. 9.5): alpha_s = 12pi/(b0 L)*(1 - b1 lnL/L), L=ln(mu^2/Lambda_nf^2),
  // b0=33-2nf, b1=2(153-19nf)/(33-2nf)^2. Lambda_{5,4,3} matched to alpha_s(M_Z)=0.1365 (mc=1.5,
  // mb=4.8) -- the exact iterated values Pythia computes (validated bit-identical, max relErr 0).
  const double L5_2=0.55304629*0.55304629, L4_2=0.72293027*0.72293027, L3_2=0.75652039*0.75652039;
  double s2min=1.33*1.33*L3_2; if(mu2<s2min) mu2=s2min;     // Pythia SAFETYMARGIN2 low-scale freeze
  double Lam2,b0,b1;
  if(mu2>mb2){ Lam2=L5_2; b0=23.0; b1=348.0/529.0; }
  else if(mu2>mc2){ Lam2=L4_2; b0=25.0; b1=462.0/625.0; }
  else { Lam2=L3_2; b0=27.0; b1=64.0/81.0; }
  double ls=log(mu2/Lam2);
  double a=12.0*M_PI/(b0*ls)*(1.0-b1*log(ls)/ls);
#else
  const double mZ2=MZ*MZ;
  const double b5=(33.0-2.0*5.0)/(12.0*M_PI), b4=(33.0-2.0*4.0)/(12.0*M_PI), b3=(33.0-2.0*3.0)/(12.0*M_PI);
  double inv=1.0/ASMZ;
  if(mu2>=mb2){ inv+=b5*log(mu2/mZ2); }
  else { inv+=b5*log(mb2/mZ2);
    if(mu2>=mc2){ inv+=b4*log(mu2/mb2); }
    else { inv+=b4*log(mc2/mb2)+b3*log(mu2/mc2); } }
  double a=1.0/inv;
#endif
#ifdef USE_CMW
  // Catani-Marchesini-Webber rescaling for soft-gluon coherence (NLL accuracy):
  // alpha_s -> alpha_s (1 + alpha_s K/2pi),  K = C_A(67/18 - pi^2/6) - 5 n_f/9.
  { int nf=(mu2>=mb2)?5:((mu2>=mc2)?4:3);
    double K=3.0*(67.0/18.0 - M_PI*M_PI/6.0) - 5.0*nf/9.0;
    a*=(1.0 + a*K/(2.0*M_PI)); }
#endif
  return (a>0.0 && a<10.0)? a : 10.0;
}
#ifdef ME_FIRST
// Generate the HARDEST emission directly from the exact O(alpha_s) gamma*/Z->q qbar g matrix
// element (POWHEG-style), via a proper pT-ORDERED Sudakov veto (a one-shot sample would miss
// the Sudakov -> wrong 2-jet/3-jet rate). Variables t=pT^2/Q^2, y=0.5 ln((1-x1)/(1-x2)):
//   dsigma ~ alpha_s(t) (x1^2+x2^2) dt/t dy,  pT^2=Q^2 (1-x1)(1-x2),  x1+x2>=1 (Dalitz).
// Overestimate (x1^2+x2^2)<=2, alpha_s(t)<=alpha_s(ptmin2); first accept = hardest w/ Sudakov.
__host__ __device__ inline bool sampleFirstEmission(uint64_t& ctr,double Q2,double ptmin2,
                                                    double* x1o,double* x2o,double* pt2o){
  double tmin=ptmin2/Q2, tmax=0.25; if(tmin>=tmax) return false;
  double ymax=acosh(1.0/(2.0*sqrt(tmin)));    // max rapidity (Dalitz boundary at tmin)
  double amax=alphaS(ptmin2);                  // SAME overestimate the shower uses
  const double CF2pi=(4.0/3.0)/(2.0*M_PI);     // C_F/2pi prefactor of the ME rate
  // Integrated overestimate per d(ln t): amax*(C_F/2pi)*2*(2 ymax). The C_F/2pi MUST be here
  // (it sets the absolute emission rate); the accept (alpha_s/amax)*(x1^2+x2^2)/2 then matches.
  double Iover=4.0*amax*ymax*CF2pi;
  double t=tmax;
  for(int it=0; it<5000; ++it){
    double R1=u01(splitmix64(ctr++)); t*=pow(R1,1.0/Iover);   // Sudakov step downward
    if(t<tmin) return false;                                  // no emission -> 2-jet
    double y =(2.0*u01(splitmix64(ctr++))-1.0)*ymax;          // draw 2: rapidity
    double Ra=u01(splitmix64(ctr++));                         // draw 3: ME-shape + running-as
    double st=sqrt(t), u=st*exp(y), v=st*exp(-y);
    if(u>=1.0||v>=1.0||u+v>1.0) continue;                     // outside the Dalitz region
    double X1=1.0-u, X2=1.0-v, p2=t*Q2;
    if(Ra < (alphaS(p2)/amax)*0.5*(X1*X1+X2*X2)){ *x1o=X1;*x2o=X2;*pt2o=p2; return true; }
  }
  return false;
}
#endif
__host__ __device__ inline double mass2(const double* a,const double* b){
  double e=a[3]+b[3],x=a[0]+b[0],y=a[1]+b[1],z=a[2]+b[2]; return e*e-x*x-y*y-z*z;
}
__host__ __device__ inline void cp4(double* d,const double* s){ d[0]=s[0];d[1]=s[1];d[2]=s[2];d[3]=s[3]; }
__host__ __device__ inline void boostBy(const double* q,double ex,double ey,double ez,double gamma,double* o){
  double bdq=ex*q[0]+ey*q[1]+ez*q[2], e2=ex*ex+ey*ey+ez*ez;
  double k=(e2>1e-18)?((gamma-1.0)*bdq/e2+gamma*q[3]):0.0;
  o[0]=q[0]+k*ex; o[1]=q[1]+k*ey; o[2]=q[2]+k*ez; o[3]=gamma*(q[3]+bdq);
}
__host__ __device__ inline void rotZto(double nx,double ny,double nz,double* v){
  double cz=(nz>1?1:(nz<-1?-1:nz)); double th=acos(cz), ph=atan2(ny,nx);
  double ct=cos(th),st=sin(th),cp=cos(ph),sp=sin(ph);
  double x=v[0],y=v[1],z=v[2];
  double x1=x*ct+z*st, y1=y, z1=-x*st+z*ct;
  v[0]=x1*cp-y1*sp; v[1]=x1*sp+y1*cp; v[2]=z1;
}
__host__ __device__ inline bool doKin(const double* R,const double* C,double pT2,double z,double phi,
                                      double* oR,double* oE,double* oC){
  double m2Dip=mass2(R,C); if(m2Dip<=0) return false; double mDip=sqrt(m2Dip);
  double m2=pT2/(z*(1.0-z));               if(m2>=m2Dip) return false;
  double eRPE =0.5*(m2Dip+m2)/mDip, pzRPE=0.5*(m2Dip-m2)/mDip;
  double pT2c =m2*(eRPE*eRPE*z*(1.0-z)-0.25*m2)/(pzRPE*pzRPE); if(pT2c<0) return false;
  double pTc=sqrt(pT2c);
  double pzRad=(eRPE*eRPE*z-0.5*m2)/pzRPE, pzEmt=(eRPE*eRPE*(1.0-z)-0.5*m2)/pzRPE;
  double cph=cos(phi),sph=sin(phi);
  double qR[4]={ pTc*cph, pTc*sph, pzRad, sqrt(pTc*pTc+pzRad*pzRad)};
  double qE[4]={-pTc*cph,-pTc*sph, pzEmt, sqrt(pTc*pTc+pzEmt*pzEmt)};
  double qC[4]={0.0,0.0,-pzRPE, pzRPE};
  double PE=R[3]+C[3], gamma=PE/mDip, bx=(R[0]+C[0])/PE,by=(R[1]+C[1])/PE,bz=(R[2]+C[2])/PE;
  double Rcm[4]; boostBy(R,-bx,-by,-bz,gamma,Rcm);
  double nn=sqrt(Rcm[0]*Rcm[0]+Rcm[1]*Rcm[1]+Rcm[2]*Rcm[2]);
  double nx,ny,nz; if(nn<1e-12){nx=0;ny=0;nz=1;}else{nx=Rcm[0]/nn;ny=Rcm[1]/nn;nz=Rcm[2]/nn;}
  rotZto(nx,ny,nz,qR); rotZto(nx,ny,nz,qE); rotZto(nx,ny,nz,qC);
  boostBy(qR,bx,by,bz,gamma,oR); boostBy(qE,bx,by,bz,gamma,oE); boostBy(qC,bx,by,bz,gamma,oC);
  return true;
}
// Shower one event into local arrays P (px,py,pz,E) and id; return nPartons.
__host__ __device__ __noinline__ int showerEvent(double* P,int* id,uint64_t ctr){  // __noinline__: keep the kernel compilable (don't inline this giant fn); physics identical
  int q0=1;
#ifdef ZFLAV
  // Initial qqbar flavour drawn from the Z hadronic branching fractions (PDG):
  // d 0.2203, u 0.1709, s 0.2203, c 0.1706, b 0.2179 (sum 1). Default build (no ZFLAV): always d.
  // NOTE: the parton shower here is massless, so b/c get no dead-cone suppression (over-radiates a
  // little at small angle); the dominant Z-flavour effect on track yield is B/D *decays* (-DHFDECAY).
  { double rf=u01(splitmix64(ctr++));
    q0 = (rf<0.2203)?1 : (rf<0.3912)?2 : (rf<0.6115)?3 : (rf<0.7821)?4 : 5; }
#endif
  P[0]=0;P[1]=0;P[2]= EBEAM;P[3]=EBEAM; id[0]= q0;
  P[4]=0;P[5]=0;P[6]=-EBEAM;P[7]=EBEAM; id[1]=-q0;
  int n=2; double pT2=0.25*MZ*MZ;
#ifdef ME_FIRST
  // Replace the first/hardest emission with the EXACT Z->qqg ME, then shower below its pT
  // (color order q, g, qbar). No (1+z^2) kernel here -> no double-count with the LL veto.
  { double x1,x2,pt2;
    if(sampleFirstEmission(ctr, MZ*MZ, PT2MIN, &x1,&x2,&pt2)){
      double E1=0.5*MZ*x1, E2=0.5*MZ*x2, E3=MZ-E1-E2;
      double c12=1.0-2.0*(x1+x2-1.0)/(x1*x2); if(c12>1.0)c12=1.0; if(c12<-1.0)c12=-1.0;
      double s12=sqrt(fmax(0.0,1.0-c12*c12));
      double phi=2.0*M_PI*u01(splitmix64(ctr++)), cph=cos(phi),sph=sin(phi);
      double qbx=E2*s12, qbz=E2*c12;               // qbar; q is along +z
      double gx=-qbx, gz=-(E1+qbz);                 // g = -(q+qbar) (q has no x)
      P[0]=0;P[1]=0;P[2]=E1;P[3]=E1; id[0]=q0;                                  // q  (+z)
      P[4]=gx*cph; P[5]=gx*sph; P[6]=gz; P[7]=E3; id[1]=21;                     // g
      P[8]=qbx*cph; P[9]=qbx*sph; P[10]=qbz; P[11]=E2; id[2]=-q0;               // qbar
      n=3; pT2=pt2;
    }
  }
#endif
  for(int step=0; step<MAXP; ++step){
    double bestT=PT2MIN, bestZ=0; int bRad=-1,bRec=-1;
#ifdef GLUON_SPLIT
    bool bSplit=false; int bFlav=0;     // g->qqbar: chosen channel + quark flavour of the best trial
#endif
    for(int i=0;i<n;++i){
      for(int side=0;side<2;++side){
        int rec; bool valid;
        if(side==0){ rec=i+1; valid=(rec<n)&&(id[i]>0||id[i]==21); }
        else       { rec=i-1; valid=(rec>=0)&&(id[i]<0||id[i]==21); }
        if(!valid) continue;
        double m2Dip=mass2(P+4*i,P+4*rec); if(m2Dip<=4.0*PT2MIN) continue;
        bool isQ=(id[i]!=21); double colFac=isQ?(4.0/3.0):(3.0/2.0);
        double zc=0.5-sqrt(fmax(0.0,0.25-PT2MIN/m2Dip)); if(zc<1e-10) zc=PT2MIN/m2Dip;
        if(zc>0.499) continue; double zmaxc=1.0-zc;
        double amax=alphaS(PT2MIN);
        double Iover=(amax/(2.0*M_PI))*colFac*2.0*log(zmaxc/zc); if(Iover<=0) continue;
#ifdef GLUON_SPLIT
        // g->qqbar trial channel (gluon radiators only): flat-z, kernel T_R[z^2+(1-z)^2], T_R=1/2.
        // nfTrial=5 (the MAX) so the flat overestimate bounds the active-flavour sum at every scale.
        // A gluon sits in TWO dipoles and each end generates the trial; unlike g->gg (genuine
        // coherent radiation off each colour line) the g->qqbar CONVERSION of the whole gluon must
        // be SHARED 1/2 per end or it double counts -> coefficient T_R*(1/2)=0.5*0.5=0.25. (Verified:
        // without the 1/2, N_gqq came out 2.14x the Pythia reference.) Adding IoverSplit only shrinks
        // the trial step; the veto stays EXACT for g->gg (unbiased).
        double IoverSplit=(!isQ)?(amax/(2.0*M_PI))*0.25*5.0*(zmaxc-zc):0.0;
#ifdef G2QQ_WEIGHT4
        IoverSplit*=G2Q4_OVER;   // oversample: the option-4 weight exceeds 1 in places -> keep the veto valid
#endif
        Iover+=IoverSplit;
#endif
        double t=pT2;
        for(int it=0; it<20000; ++it){
          double R1=u01(splitmix64(ctr++)); t*=pow(R1,1.0/Iover);
          if(t<PT2MIN) break;
#ifdef GLUON_SPLIT
          // FIXED 5 draws/iteration for EVERY dipole end (gluon or quark): Rc,Rf are drawn even for
          // quark radiators (where they go unused) so the per-thread RNG phase never depends on
          // parton flavour -> control flow stays 100% GPU==CPU. (Off: 3 draws, byte-identical.)
          double Rc=u01(splitmix64(ctr++));                              // channel pick
          double Rf=u01(splitmix64(ctr++));                             // flavour pick
          double R2=u01(splitmix64(ctr++));
          bool doSplit=(!isQ)&&(Rc*Iover<IoverSplit);
          double zz=doSplit ? (zc+(zmaxc-zc)*R2)                        // flat-z (kernel is pole-free)
                            : (1.0-(1.0-zc)*pow(zc/(1.0-zc),R2));       // log-z (emit, unchanged)
#else
          double R2=u01(splitmix64(ctr++));
          double zz=1.0-(1.0-zc)*pow(zc/(1.0-zc),R2);
#endif
          double R3=u01(splitmix64(ctr++));
          double zmn=0.5-sqrt(fmax(0.0,0.25-t/m2Dip)), zmx=1.0-zmn;
          double m2v=t/(zz*(1.0-zz));
          if(zz>zmn && zz<zmx && m2v<m2Dip){
            double dal=zz*(1.0-zz)*(m2Dip+m2v)*(m2Dip+m2v);
            if(m2v*m2Dip<dal){
#ifdef GLUON_SPLIT
              double w; int flav=0;
              if(doSplit){
                flav=1+(int)(5.0*Rf); if(flav>5)flav=5;                // uniform flavour {1..nGluonToQuark=5}
                // Pythia quark m0: u,d=0.33, s=0.50, c=1.50, b=4.80 (ParticleData.xml). The pair
                // threshold m2_qq >= THRESHM2*m_q^2 (4.004) and betaQ both suppress low-mass pairs;
                // treating u,d as massless overproduced uds ~9% -> use m0 like Pythia. No nf-count
                // veto (alphaS() handles the running nf separately); the kinematic threshold IS the gate.
                double mq=(flav==5)?4.8:((flav==4)?1.5:((flav==3)?0.5:0.33));
                if(m2v>=4.004*mq*mq){                                  // = Pythia THRESHM2 pair threshold
                  double beta=sqrt(fmax(0.0,1.0-4.0*mq*mq/m2v));       // = Pythia betaQ
#ifdef G2QQ_WEIGHT4
                  // Pythia DEFAULT weightGluonToQuark=4 (SimpleTimeShower.cc 2937-2943): the zCosThe
                  // reshape + (1+m2Rat)/(1-m2Rat)*pow3(1-m2Rat) = *(1+m2Rat)(1-m2Rat)^2 damping.
                  // zCosThe pushes z outside [0,1] so w>1 in places -> the IoverSplit*G2Q4_OVER
                  // oversample (and /G2Q4_OVER here) keep accept<=1, rate unchanged. Matches option-4's
                  // RATE; produced kinematics stay massless (the zCosThe massive-recoil construction is
                  // future) -> a rate-level match to Pythia's DEFAULT, not its exact z-distribution.
                  double m2Rat=m2v/m2Dip, ratioQ=mq*mq/m2v;
                  double zct=((1.0+m2Rat)*zz-m2Rat)/(1.0-m2Rat);
                  w=beta*(zct*zct+(1.0-zct)*(1.0-zct)+8.0*ratioQ*zct*(1.0-zct))
                      *(1.0+m2Rat)*(1.0-m2Rat)*(1.0-m2Rat)/G2Q4_OVER;
#else
                  w=beta*(zz*zz+(1.0-zz)*(1.0-zz));                    // option 1: plain DGLAP (NO T_R; in IoverSplit)
#endif
                } else w=-1.0;                                         // below pair threshold -> veto
              } else w=isQ?0.5*(1.0+zz*zz):0.5*(1.0+zz*zz*zz);
              if(w>0.0 && R3<(alphaS(t)/amax)*w){
                if(t>bestT){ bestT=t; bestZ=zz; bRad=i; bRec=rec; bSplit=doSplit; bFlav=flav; } break; }
#else
              double w=isQ?0.5*(1.0+zz*zz):0.5*(1.0+zz*zz*zz);
              if(R3<(alphaS(t)/amax)*w){ if(t>bestT){ bestT=t; bestZ=zz; bRad=i; bRec=rec; } break; }
#endif
            }
          }
        }
      }
    }
    if(bRad<0) break;
    pT2=bestT;
    double phi=2.0*M_PI*u01(splitmix64(ctr++));
    double oR[4],oE[4],oC[4];
    if(!doKin(P+4*bRad,P+4*bRec,bestT,bestZ,phi,oR,oE,oC)) continue;
    if(n>=MAXP) break;
    if(bRec==bRad+1){
      for(int k=n;k>bRad+1;--k){ cp4(P+4*k,P+4*(k-1)); id[k]=id[k-1]; }
      cp4(P+4*bRad,oR); cp4(P+4*(bRad+1),oE); cp4(P+4*(bRad+2),oC);
#ifdef GLUON_SPLIT
      if(bSplit){ id[bRad]=-bFlav; id[bRad+1]=bFlav; }  // qbar' (connects left), q' (connects right)
      else id[bRad+1]=21;
#else
      id[bRad+1]=21;
#endif
    } else {
      for(int k=n;k>bRad;--k){ cp4(P+4*k,P+4*(k-1)); id[k]=id[k-1]; }
      cp4(P+4*(bRad-1),oC); cp4(P+4*bRad,oE); cp4(P+4*(bRad+1),oR);
#ifdef GLUON_SPLIT
      if(bSplit){ id[bRad]=-bFlav; id[bRad+1]=bFlav; }  // qbar' (connects left), q' (connects right)
      else id[bRad]=21;
#else
      id[bRad]=21;
#endif
    }
    n++;
  }
  return n;
}
