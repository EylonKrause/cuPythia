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

// 1-loop running alpha_s, flavour-threshold matched (n_f=5,4,3 across m_b,m_c).
__host__ __device__ inline double alphaS(double mu2){
  const double mc2=1.5*1.5, mb2=4.8*4.8, mZ2=MZ*MZ;
  const double b5=(33.0-2.0*5.0)/(12.0*M_PI), b4=(33.0-2.0*4.0)/(12.0*M_PI), b3=(33.0-2.0*3.0)/(12.0*M_PI);
  double inv=1.0/ASMZ;
  if(mu2>=mb2){ inv+=b5*log(mu2/mZ2); }
  else { inv+=b5*log(mb2/mZ2);
    if(mu2>=mc2){ inv+=b4*log(mu2/mb2); }
    else { inv+=b4*log(mc2/mb2)+b3*log(mu2/mc2); } }
  double a=1.0/inv;
  return (a>0.0 && a<10.0)? a : 10.0;
}
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
__host__ __device__ inline int showerEvent(double* P,int* id,uint64_t ctr){
  P[0]=0;P[1]=0;P[2]= EBEAM;P[3]=EBEAM; id[0]= 1;
  P[4]=0;P[5]=0;P[6]=-EBEAM;P[7]=EBEAM; id[1]=-1;
  int n=2; double pT2=0.25*MZ*MZ;
  for(int step=0; step<MAXP; ++step){
    double bestT=PT2MIN, bestZ=0; int bRad=-1,bRec=-1;
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
        double t=pT2;
        for(int it=0; it<20000; ++it){
          double R1=u01(splitmix64(ctr++)); t*=pow(R1,1.0/Iover);
          if(t<PT2MIN) break;
          double R2=u01(splitmix64(ctr++));
          double zz=1.0-(1.0-zc)*pow(zc/(1.0-zc),R2);
          double R3=u01(splitmix64(ctr++));
          double zmn=0.5-sqrt(fmax(0.0,0.25-t/m2Dip)), zmx=1.0-zmn;
          double m2v=t/(zz*(1.0-zz));
          if(zz>zmn && zz<zmx && m2v<m2Dip){
            double dal=zz*(1.0-zz)*(m2Dip+m2v)*(m2Dip+m2v);
            if(m2v*m2Dip<dal){
              double w=isQ?0.5*(1.0+zz*zz):0.5*(1.0+zz*zz*zz);
              if(R3<(alphaS(t)/amax)*w){ if(t>bestT){ bestT=t; bestZ=zz; bRad=i; bRec=rec; } break; }
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
      cp4(P+4*bRad,oR); cp4(P+4*(bRad+1),oE); id[bRad+1]=21; cp4(P+4*(bRad+2),oC);
    } else {
      for(int k=n;k>bRad;--k){ cp4(P+4*k,P+4*(k-1)); id[k]=id[k-1]; }
      cp4(P+4*(bRad-1),oC); cp4(P+4*bRad,oE); id[bRad]=21; cp4(P+4*(bRad+1),oR);
    }
    n++;
  }
  return n;
}
