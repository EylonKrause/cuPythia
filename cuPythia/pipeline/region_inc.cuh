// cuPythia hadronization — multi-region (gluon-kinked string) FOUNDATION.
// Faithful port of Pythia 8.317 StringRegion::setUp / project / pHad and StringSystem::setUp
// (FragmentationSystems.cc:535-678, .h:167-211): the light-cone region table for a gluon-
// kinked string. A gluon's momentum is shared 50/50 between the two regions it borders
// (the kink). This header is the de-risked base for the multi-region kinematicsHadron port;
// region_test.cu validates the construction's math (lightlike basis, orthonormal transverse
// axes, projection round-trip) standalone before it is wired into the fragmentation chain.
#pragma once
#include <cmath>
#include <cstdint>

#define MAXPART 64
#define MAXREG  ((MAXPART-1)*MAXPART/2)
#define R_MJOIN 0.1
#define R_TINY  1e-20

// Vec4 = {px,py,pz,E}; Minkowski dot (+ - - -).
__host__ __device__ inline double v4dot(const double* a,const double* b){
  return a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2];
}
struct Region { double pPos[4],pNeg[4],eX[4],eY[4],w2; int colPos,colNeg; bool isSetUp,isEmpty; };

// Port of StringRegion::setUp. isMassless: incoming guaranteed lightlike (cross regions).
__host__ __device__ __noinline__ void regionSetUp(Region& R,const double* p1in,const double* p2in,
                                            int col1,int col2,bool isMassless){
  const double m2Join=R_MJOIN*R_MJOIN;
  double p1[4],p2[4]; for(int k=0;k<4;++k){p1[k]=p1in[k];p2[k]=p2in[k];}
  R.isSetUp=false; R.isEmpty=false;
  // Zero the basis so an EMPTY region's pHad gives 0 (Pythia's default-constructed Vec4),
  // not uninitialized garbage — essential for determinism in the multi-region stepping.
  for(int k=0;k<4;++k){R.pPos[k]=0.0;R.pNeg[k]=0.0;R.eX[k]=0.0;R.eY[k]=0.0;} R.w2=0.0;
  if(isMassless){
    R.w2=2.0*v4dot(p1,p2);
    if(R.w2<m2Join){R.isSetUp=true;R.isEmpty=true;return;}
    for(int k=0;k<4;++k){R.pPos[k]=p1[k];R.pNeg[k]=p2[k];}
  } else {
    double m1Sq=v4dot(p1,p1), m2Sq=v4dot(p2,p2), p1p2=v4dot(p1,p2);
    R.w2=m1Sq+2.0*p1p2+m2Sq; double rootSq=p1p2*p1p2-m1Sq*m2Sq;
    if(R.w2<=0.0||rootSq<=0.0){
      if(m1Sq<0.0)m1Sq=0.0; p1[3]=sqrt(m1Sq+p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
      if(m2Sq<0.0)m2Sq=0.0; p2[3]=sqrt(m2Sq+p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
      p1p2=v4dot(p1,p2); R.w2=m1Sq+2.0*p1p2+m2Sq; rootSq=p1p2*p1p2-m1Sq*m2Sq;
    }
    if(R.w2<m2Join){R.isSetUp=true;R.isEmpty=true;return;}
    double root=sqrt(fmax(R_TINY,rootSq));
    double k1=0.5*((m2Sq+p1p2)/root-1.0), k2=0.5*((m1Sq+p1p2)/root-1.0);
    for(int k=0;k<4;++k){ R.pPos[k]=(1.0+k1)*p1[k]-k2*p2[k]; R.pNeg[k]=(1.0+k2)*p2[k]-k1*p1[k]; }
    if(R.pPos[3]<R_TINY||R.pNeg[3]<R_TINY){R.isSetUp=true;R.isEmpty=true;return;}
  }
  // Transverse axes: trial directions then Gram-Schmidt (StringRegion::setUp:588-619).
  double ePos[4],eNeg[4]; for(int k=0;k<4;++k){ePos[k]=R.pPos[k]/R.pPos[3]; eNeg[k]=R.pNeg[k]/R.pNeg[3];}
  double eDx=(ePos[0]-eNeg[0])*(ePos[0]-eNeg[0]);
  double eDy=(ePos[1]-eNeg[1])*(ePos[1]-eNeg[1]);
  double eDz=(ePos[2]-eNeg[2])*(ePos[2]-eNeg[2]);
  double eX[4]={0,0,0,0}, eY[4]={0,0,0,0};
  if(eDx<fmin(eDy,eDz)){ eX[0]=1.0; if(eDy<eDz)eY[1]=1.0; else eY[2]=1.0; }
  else if(eDy<eDz){ eX[1]=1.0; if(eDx<eDz)eY[0]=1.0; else eY[2]=1.0; }
  else { eX[2]=1.0; if(eDx<eDy)eY[0]=1.0; else eY[1]=1.0; }
  double pPosNeg=v4dot(R.pPos,R.pNeg);
  double kXPos=v4dot(eX,R.pPos)/pPosNeg, kXNeg=v4dot(eX,R.pNeg)/pPosNeg;
  double kXtmp=1.0+2.0*kXPos*kXNeg*pPosNeg; if(kXtmp<R_TINY){R.isSetUp=true;R.isEmpty=true;return;}
  double kXX=1.0/sqrt(kXtmp);
  double kYPos=v4dot(eY,R.pPos)/pPosNeg, kYNeg=v4dot(eY,R.pNeg)/pPosNeg;
  double kYX=kXX*(kXPos*kYNeg+kXNeg*kYPos)*pPosNeg;
  double kYtmp=1.0+2.0*kYPos*kYNeg*pPosNeg-kYX*kYX; if(kYtmp<R_TINY){R.isSetUp=true;R.isEmpty=true;return;}
  double kYY=1.0/sqrt(kYtmp);
  for(int k=0;k<4;++k) R.eX[k]=kXX*(eX[k]-kXNeg*R.pPos[k]-kXPos*R.pNeg[k]);
  for(int k=0;k<4;++k) R.eY[k]=kYY*(eY[k]-kYNeg*R.pPos[k]-kYPos*R.pNeg[k]-kYX*R.eX[k]);
  R.colPos=col1; R.colNeg=col2; R.isSetUp=true; R.isEmpty=false;
}
// pHad and project (StringRegion::pHad/project).
__host__ __device__ inline void regionPHad(const Region& R,double xP,double xN,double px,double py,double* o){
  for(int k=0;k<4;++k) o[k]=xP*R.pPos[k]+xN*R.pNeg[k]+px*R.eX[k]+py*R.eY[k];
}
__host__ __device__ inline void regionProject(const Region& R,const double* pIn,
                                              double& xP,double& xN,double& px,double& py){
  xP=2.0*v4dot(pIn,R.pNeg)/R.w2; xN=2.0*v4dot(pIn,R.pPos)/R.w2;
  px=-v4dot(pIn,R.eX); py=-v4dot(pIn,R.eY);
}

struct StringSys { Region reg[MAXREG]; int sizeStr,iMax,indxReg,nReg; };
__host__ __device__ inline int sysIReg(const StringSys& S,int iPos,int iNeg){
  return (iPos*(S.indxReg-iPos))/2+iNeg;
}
// Build the lowest-lying regions from the colour-ordered parton chain P (px,py,pz,E per
// parton), id (21=gluon), with sequential colour tags. Cross regions are built lazily.
__host__ __device__ inline void sysSetUp(StringSys& S,const double* P,const int* id,int n){
  S.sizeStr=n-1; S.iMax=S.sizeStr-1; S.indxReg=2*S.sizeStr+1; S.nReg=S.sizeStr*(S.sizeStr+1)/2;
  for(int r=0;r<S.nReg;++r) S.reg[r].isSetUp=false;
  for(int i=0;i<S.sizeStr;++i){
    double p1[4],p2[4];
    for(int k=0;k<4;++k){ p1[k]=P[4*i+k]*((id[i]==21)?0.5:1.0); p2[k]=P[4*(i+1)+k]*((id[i+1]==21)?0.5:1.0); }
    int col=101+i;   // forward chain (q endpoint has colour); colPos/colNeg bookkeeping only
    regionSetUp(S.reg[sysIReg(S,i,S.iMax-i)], p1,p2, col,col, false);
  }
}
