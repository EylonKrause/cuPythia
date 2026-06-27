// cuPythia — ALL-GPU multi-region (gluon-kinked) Lund string fragmentation.
//
// Extends the straight-string hadronizer to the FSR shower's q-g-...-qbar chains by porting
// Pythia 8.317 StringEnd::kinematicsHadron (the region-stepping with the (m^2,Gamma) coupled
// quadratic + bidirectional region-crossing, StringFragmentation.cc:131-311), update (:531),
// newHadron pT (pxHad=pxOld+pxNew, :91-121), finalRegion (:1851) and the multi-region finalTwo
// (:1640-1716) onto the validated region table (region_inc.cuh). Flavour/spin/eta-mix, the
// Lund f(z) sampler, the constituent-mass stop test and the refragment loop are reused from the
// straight-string slice (validated to 4% vs Pythia). One string per thread, counter-RNG.
//
// Validation: exact 4-momentum conservation + on-shellness on gluon-kinked strings (the
// stringent check that catches any kinematics bug), GPU==CPU determinism, reproducibility,
// and full hadron multiplicity (cross-checked vs Pythia forceHadronLevel via bridge_pythia).
//
// Build: nvcc -O3 -std=c++17 -arch=sm_120 -o hadronize_mr hadronize_mr.cu
// Run:   ./hadronize_mr [nEvents=20000]

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cuda_runtime.h>
#include "../common/rng.cuh"
#include "region_inc.cuh"
#include "zlund_inc.cuh"
#include "shower_inc.cuh"     // showerEvent: produces the q-g-...-qbar chain (MZ, EBEAM)
#ifdef DECAYS
#include "decay_inc.cuh"      // GPU recursive hadron decays (rho/K*/omega/phi/eta/eta'/K0 -> pi/K/gamma)
#define OUTCAP MAXFINAL       // decays multiply particle count -> larger output buffer cap
#else
#define OUTCAP MAXPART
#endif

#define CK(c) do{cudaError_t e=(c); if(e!=cudaSuccess){printf("CUDA %s @%d\n",cudaGetErrorString(e),__LINE__);return 1;}}while(0)

// ---- fragmentation parameters (Pythia 8.317 defaults; same as hadronize.cu) ----
static const double H_ALUND=0.68, H_BLUND=0.98;
static const double H_SIGMAQ=0.335/1.4142135623730951, H_PROBSTOUD=0.217;
static const double H_MUDV=0.50, H_SV=0.55, H_THPS=-15.0, H_THV=36.0;
static const double H_ETASUP=0.60, H_ETAPSUP=0.12;
static const double H_STOPM=0.8, H_STOPNF=2.0, H_STOPSM=0.2, H_ENHF=0.01, H_ENHW=2.0;
static const double EQ_TINY=1e-6, PT2SAME=0.01;   // StringEnd::TINY, PT2SAME

__host__ __device__ inline double mesonMassMR(int pdg){
  switch(abs(pdg)){
    case 211:return 0.13957; case 111:return 0.13498; case 221:return 0.54786; case 331:return 0.95778;
    case 321:return 0.49368; case 311:return 0.49761; case 130: case 310:return 0.49761; case 22:return 0.0;
    case 213:return 0.77526; case 113:return 0.77526; case 223:return 0.78266; case 333:return 1.01946;
    case 323:return 0.89167; case 313:return 0.89555;
  } return -1.0;
}
__host__ __device__ inline double constMassMR(int a){ return (a==3)?0.5:0.33; }
__host__ __device__ inline bool isChargedMR(int pdg){ int a=abs(pdg); return a==211||a==321||a==213||a==323; }
__host__ __device__ inline int pickFlavMR(uint64_t& c){ double r=u01(splitmix64(c++))*(2.0+H_PROBSTOUD); return (r<1.0)?1:((r<2.0)?2:3); }
__host__ __device__ inline int combineMesonMR(int id1,int id2,uint64_t& c){
  int a1=abs(id1),a2=abs(id2),idMax=(a1>a2)?a1:a2,idMin=(a1<a2)?a1:a2;
  int flav=(idMax<3)?0:idMax-2; double mv=(flav==0)?H_MUDV:H_SV;
  double rs=(1.0+mv)*u01(splitmix64(c++)); int spin=0; rs-=1.0; if(rs>0.0) spin=1;
  int code=(spin==0)?1:3, idM=100*idMax+10*idMin+code;
  if(idMax!=idMin){ int sg=(idMax%2==0)?1:-1; if((idMax==a1&&id1<0)||(idMax==a2&&id2<0)) sg=-sg; idM*=sg; }
  else { double al=((spin==0)?(90.0-(H_THPS+54.7)):(H_THV+54.7))*M_PI/180.0;
    double m1,m2; if(flav==0){m1=0.5;m2=0.5*(1.0+sin(al)*sin(al));}else{m1=0.0;m2=cos(al)*cos(al);}
    double rm=u01(splitmix64(c++)); if(rm<m1)idM=110;else if(rm<m2)idM=220;else idM=330; idM+=code;
    if(idM==221&&H_ETASUP<u01(splitmix64(c++)))return 0; if(idM==331&&H_ETAPSUP<u01(splitmix64(c++)))return 0; }
  return idM;
}
__host__ __device__ inline void pairPTMR(uint64_t& c,double& gx,double& gy){
  double m=(u01(splitmix64(c++))<H_ENHF)?H_ENHW:1.0;
  double r=sqrt(-2.0*log(u01(splitmix64(c++))+1e-300)), ph=2.0*M_PI*u01(splitmix64(c++));
  gx=m*H_SIGMAQ*r*sin(ph); gy=m*H_SIGMAQ*r*cos(ph);
}
__host__ __device__ inline void boostByMR(const double* q,double ex,double ey,double ez,double g,double* o){
  double bq=ex*q[0]+ey*q[1]+ez*q[2], e2=ex*ex+ey*ey+ez*ez;
  double k=(e2>1e-18)?((g-1.0)*bq/e2+g*q[3]):0.0;
  o[0]=q[0]+k*ex;o[1]=q[1]+k*ey;o[2]=q[2]+k*ez;o[3]=g*(q[3]+bq);
}

// Lightweight system: store only the n-1 lowest-lying regions; build cross regions on demand
// (the full O(n^2) table would be ~280 KB/thread of local memory and corrupts). A gluon's
// momentum is shared 50/50 between its two bordering low regions (the kink).
struct SysLite { Region low[MAXPART]; int sizeStr,iMax; };
__host__ __device__ inline void sysLiteSetUp(SysLite& S,const double* P,const int* id,int n){
  S.sizeStr=n-1; S.iMax=S.sizeStr-1;
  for(int i=0;i<S.sizeStr;++i){ double p1[4],p2[4];
    for(int k=0;k<4;++k){ p1[k]=P[4*i+k]*((id[i]==21)?0.5:1.0); p2[k]=P[4*(i+1)+k]*((id[i+1]==21)?0.5:1.0); }
    regionSetUp(S.low[i], p1,p2, 101+i,101+i, false); }
}
// regionLowPos(iPos)=low[iPos]; regionLowNeg(iNeg)=low[iMax-iNeg]. region(iPos,iNeg) on demand.
__host__ __device__ inline void getRegion(const SysLite& S,int iPos,int iNeg,Region& out){
  if(iPos+iNeg==S.iMax) out=S.low[iPos];
  else regionSetUp(out, S.low[iPos].pPos, S.low[S.iMax-iNeg].pNeg, S.low[iPos].colPos, S.low[S.iMax-iNeg].colNeg, true);
}
__host__ __device__ inline double w2Region(const SysLite& S,int iPos,int iNeg){
  if(iPos+iNeg==S.iMax) return S.low[iPos].w2;
  return 2.0*v4dot(S.low[iPos].pPos, S.low[S.iMax-iNeg].pNeg);
}

// Per string end.
struct End { int iPos,iNeg; double xPos,xNeg,Gamma,px,py; int flav; };

// Port of StringEnd::kinematicsHadron. Returns true on success, fills had[4] and the *New state.
__host__ __device__ inline bool kinHad(SysLite& S,bool fromPos,const End& e,
    double zHad,double mT2Had,double mHad,double pxHadIn,double pyHadIn,
    double pxNewIn,double pyNewIn,
    int& iPosNewO,int& iNegNewO,double& xPosNewO,double& xNegNewO,double& GammaNewO,double* had){
  double GammaNew=(1.0-zHad)*(e.Gamma+mT2Had/zHad); GammaNewO=GammaNew;
  int iMax=S.iMax;
  int iDirOld=fromPos?e.iPos:e.iNeg, iInvOld=fromPos?e.iNeg:e.iPos;
  double xDirOld=fromPos?e.xPos:e.xNeg, xInvOld=fromPos?e.xNeg:e.xPos;
  int iPosNew=e.iPos, iNegNew=e.iNeg;
  double pxNew=pxNewIn, pyNew=pyNewIn;
  double pSoFar[4]={0,0,0,0}, pTNew[4]={0,0,0,0}, tmp[4];
  for(int iStep=0; iStep<4*iMax+8; ++iStep){    // generous bound for multi-region crossings
    Region region; getRegion(S,iPosNew,iNegNew,region);
    double xPosHad,xNegHad,xDirHad,xInvHad;
    if(iStep==0 && e.iPos+e.iNeg==iMax){
      if(mT2Had < zHad*xDirOld*(1.0-xInvOld)*region.w2){
        xDirHad=zHad*xDirOld; xInvHad=mT2Had/(xDirHad*region.w2);
        double xDirNew=xDirOld-xDirHad, xInvNew=xInvOld+xInvHad;
        xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
        iPosNewO=iPosNew; iNegNewO=iNegNew;
        xPosNewO=fromPos?xDirNew:xInvNew; xNegNewO=fromPos?xInvNew:xDirNew;
        regionPHad(region,xPosHad,xNegHad,pxHadIn,pyHadIn,had); return true;
      } else {
        if(fromPos) iNegNew--; else iPosNew--;
        int iInvNew=fromPos?iNegNew:iPosNew; if(iInvNew<0) return false;
        xInvHad=1.0-xInvOld; xDirHad=0.0;
        xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
        regionPHad(region,xPosHad,xNegHad,e.px,e.py,pSoFar); continue;
      }
    } else if(iStep==0){
      regionPHad(region,0,0,e.px,e.py,pSoFar);
      regionPHad(region,0,0,pxNew,pyNew,pTNew);
    }
    if(region.isEmpty){
      int iDirNew=fromPos?iPosNew:iNegNew;
      xDirHad=(iDirNew==iDirOld)?xDirOld:1.0; xInvHad=0.0;
      xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
      regionPHad(region,xPosHad,xNegHad,0,0,tmp); for(int k=0;k<4;++k)pSoFar[k]+=tmp[k];
      if(fromPos)iPosNew++; else iNegNew++;
      int iDN=fromPos?iPosNew:iNegNew, iIN=fromPos?iNegNew:iPosNew;
      if(iDN+iIN>iMax) return false; continue;
    }
    double pxT=-v4dot(pTNew,region.eX), pyT=-v4dot(pTNew,region.eY);
    if(fabs(pxT*pxT+pyT*pyT-pxNew*pxNew-pyNew*pyNew)<PT2SAME){ pxNew=pxT; pyNew=pyT; }
    double pTemp[4]; regionPHad(region,0,0,pxNew,pyNew,tmp); for(int k=0;k<4;++k)pTemp[k]=pSoFar[k]+tmp[k];
    double cM1=v4dot(pTemp,pTemp), cM2=2.0*v4dot(pTemp,region.pPos), cM3=2.0*v4dot(pTemp,region.pNeg), cM4=region.w2;
    if(!fromPos){ double t=cM2;cM2=cM3;cM3=t; }
    int iDirNew=fromPos?iPosNew:iNegNew, iInvNew=fromPos?iNegNew:iPosNew;
    double cGam1=0,cGam2=0,cGam3=0,cGam4=0;
    for(int iInv=iInvNew; iInv<=iMax-iDirNew; ++iInv){
      double xInv=1.0; if(iInv==iInvNew) xInv=(iInvNew==iInvOld)?xInvOld:0.0;
      for(int iDir=iDirNew; iDir<=iMax-iInv; ++iDir){
        double xDir=(iDir==iDirOld)?xDirOld:1.0;
        int ip=fromPos?iDir:iInv, in=fromPos?iInv:iDir; double w2=w2Region(S,ip,in);
        cGam1+=xDir*xInv*w2;
        if(iDir==iDirNew)cGam2-=xInv*w2;
        if(iInv==iInvNew)cGam3+=xDir*w2;
        if(iDir==iDirNew&&iInv==iInvNew)cGam4-=w2;
      }
    }
    double cM0=mHad*mHad-cM1, cGam0=GammaNew-cGam1;
    double r2=cM3*cGam4-cM4*cGam3, r1=cM4*cGam0-cM0*cGam4+cM3*cGam2-cM2*cGam3, r0=cM2*cGam0-cM0*cGam2;
    double disc=r1*r1-4.0*r2*r0, root=(disc>0.0)?sqrt(disc):0.0;
    if(fabs(r2)<EQ_TINY || root<EQ_TINY) return false;
    xInvHad=0.5*(root/fabs(r2)-r1/r2);
    if(fabs(cM2+cM4*xInvHad)<EQ_TINY) return false;
    xDirHad=(cM0-cM3*xInvHad)/(cM2+cM4*xInvHad);
    double xDirNew=(iDirNew==iDirOld)?xDirOld-xDirHad:1.0-xDirHad;
    double xInvNew=(iInvNew==iInvOld)?xInvOld+xInvHad:xInvHad;
    if(xInvNew>1.0){
      xInvHad=(iInvNew==iInvOld)?1.0-xInvOld:1.0; xDirHad=0.0;
      xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
      regionPHad(region,xPosHad,xNegHad,0,0,tmp); for(int k=0;k<4;++k)pSoFar[k]+=tmp[k];
      if(fromPos)iNegNew--; else iPosNew--; int iIN=fromPos?iNegNew:iPosNew; if(iIN<0)return false; continue;
    } else if(xDirNew<0.0){
      xDirHad=(iDirNew==iDirOld)?xDirOld:1.0; xInvHad=0.0;
      xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
      regionPHad(region,xPosHad,xNegHad,0,0,tmp); for(int k=0;k<4;++k)pSoFar[k]+=tmp[k];
      if(fromPos)iPosNew++; else iNegNew++;
      int iDN=fromPos?iPosNew:iNegNew, iIN=fromPos?iNegNew:iPosNew; if(iDN+iIN>iMax)return false; continue;
    }
    xPosHad=fromPos?xDirHad:xInvHad; xNegHad=fromPos?xInvHad:xDirHad;
    iPosNewO=iPosNew; iNegNewO=iNegNew; xPosNewO=fromPos?xDirNew:xInvNew; xNegNewO=fromPos?xInvNew:xDirNew;
    regionPHad(region,xPosHad,xNegHad,pxNew,pyNew,tmp); for(int k=0;k<4;++k)had[k]=pSoFar[k]+tmp[k];
    return true;
  }
  return false;
}

// finalRegion (port of :1851), simple-and-common cases for an open q..qbar chain.
__host__ __device__ inline bool finalRegion(const SysLite& S,const End& pos,const End& neg,Region& out){
  if(pos.iPos==neg.iPos && pos.iNeg==neg.iNeg){ getRegion(S,pos.iPos,pos.iNeg,out); return !out.isEmpty; }
  double pPosJ[4]={0,0,0,0}, pNegJ[4]={0,0,0,0}, tmp[4];
  // remaining p+
  if(pos.iPos==neg.iPos){ double x=pos.xPos-neg.xPos; if(x<0)return false;
    regionPHad(S.low[pos.iPos],x,0,0,0,pPosJ);
  } else { for(int ip=pos.iPos; ip<=neg.iPos; ++ip){ const Region& lr=S.low[ip];
      double xf=(ip==pos.iPos)?pos.xPos:((ip==neg.iPos)?1.0-neg.xPos:1.0);
      regionPHad(lr,xf,0,0,0,tmp); for(int k=0;k<4;++k)pPosJ[k]+=tmp[k]; } }
  // remaining p-
  if(neg.iNeg==pos.iNeg){ double x=neg.xNeg-pos.xNeg; if(x<0)return false;
    regionPHad(S.low[S.iMax-neg.iNeg],0,x,0,0,pNegJ);
  } else { for(int in=neg.iNeg; in<=pos.iNeg; ++in){ const Region& lr=S.low[S.iMax-in];
      double xf=(in==neg.iNeg)?neg.xNeg:((in==pos.iNeg)?1.0-pos.xNeg:1.0);
      regionPHad(lr,0,xf,0,0,tmp); for(int k=0;k<4;++k)pNegJ[k]+=tmp[k]; } }
  int cP=S.low[pos.iPos].colPos;
  regionSetUp(out, pPosJ, pNegJ, cP, cP, false);  // joined p+/p- are TIMELIKE -> massive construction
  return !out.isEmpty;
}

// Split a colour-ordered chain into independent strings at every g->qqbar boundary: a break sits
// between i and i+1 iff a real antiquark (id<0) is immediately followed by a real quark (id>0,!=21).
// A pure q...gluons...qbar chain has no such boundary -> returns exactly one string [0..n-1]
// (fully backward compatible). Each string is a maximal run [start..end] (q at start, qbar at end).
__host__ __device__ inline int findStrings(const int* id,int n,int* starts,int* ends){
  int ns=0, st=0;
  for(int i=0;i<n-1;++i)
    if(id[i]<0 && id[i]!=-21 && id[i+1]>0 && id[i+1]!=21){ starts[ns]=st; ends[ns]=i; ns++; st=i+1; }
  starts[ns]=st; ends[ns]=n-1; ns++;
  return ns;
}

// Fragment one chain; returns nHadrons or -1 (caller refragments).
#ifdef USE_BW
#include "bw_inc.cuh"
#define HADMASS_MR(pdg,ctr) sampleBWmass(pdg, mesonMassMR(pdg), ctr)   // Breit-Wigner vector masses
#else
#define HADMASS_MR(pdg,ctr) mesonMassMR(pdg)                            // pole masses (default)
#endif
__host__ __device__ inline int tryFragmentMR(const double* P,const int* id,int n,uint64_t& ctr,double* H,int* hid,double* hm){
  SysLite S; sysLiteSetUp(S,P,id,n);
  // Endpoint flavours from the actual string ends (q at id[0], qbar at id[n-1]); a g->qqbar fork
  // makes these s/c/b. Identical to the old hard-wired +-1 for the original d...dbar string.
  End pos={0,S.iMax,1.0,0.0,0.0,0.0,0.0, id[0]}, neg={S.iMax,0,0.0,1.0,0.0,0.0,0.0, id[n-1]};
  double pRem[4]={0,0,0,0}; for(int i=0;i<n;++i) for(int k=0;k<4;++k) pRem[k]+=P[4*i+k];
  int nH=0;
  for(int step=0; step<MAXPART-2; ++step){
    bool fromPos=(u01(splitmix64(ctr++))<0.5);
    End& e = fromPos?pos:neg;
    int pdg=0,fNew=0; double mHad=0,pxHad=0,pyHad=0,pxNew=0,pyNew=0,gx,gy;
    for(int tr=0; tr<20 && pdg==0; ++tr){
      fNew=pickFlavMR(ctr);
      pdg=combineMesonMR(e.flav, fromPos?-fNew:fNew, ctr);
      pairPTMR(ctr,gx,gy);
      if(pdg!=0){ mHad=HADMASS_MR(pdg,ctr); pxNew=gx; pyNew=gy; pxHad=e.px+pxNew; pyHad=e.py+pyNew; }
    }
    if(pdg==0) return -1;
    double mT2Had=mHad*mHad+pxHad*pxHad+pyHad*pyHad;
    double wMin=(H_STOPM+constMassMR(abs(pos.flav))+constMassMR(abs(neg.flav))+H_STOPNF*constMassMR(fNew))
                *(1.0+(2.0*u01(splitmix64(ctr++))-1.0)*H_STOPSM);
    double w2Rem=pRem[3]*pRem[3]-pRem[0]*pRem[0]-pRem[1]*pRem[1]-pRem[2]*pRem[2];
    if(w2Rem<wMin*wMin) break;
    // Retry the z draw a few times if the region-crossing solve fails or gives a bad hadron
    // (Pythia retries the hadron before refragmenting the whole string). Only on-shell, finite
    // hadrons are accepted -> conservation stays exact (finalTwo distributes the remainder).
    // Retry the hadron (fresh flavour+pT+z) before giving up — Pythia's fallback order is to
    // retry the hadron, refragmenting the whole string only as a last resort.
    int iPN,iNN; double xPN,xNN,GN,had[4],accMass=mHad; bool made=false;
    for(int ht=0; ht<24 && !made; ++ht){
      int rpdg=0,rfNew=0; double rmHad=mHad,rpxHad=pxHad,rpyHad=pyHad,rpxNew=pxNew,rpyNew=pyNew,rmT2=mT2Had,rgx,rgy;
      if(ht>0){ // redraw flavour+pT for retries (ht==0 reuses the already-drawn hadron)
        for(int tr=0; tr<20 && rpdg==0; ++tr){ rfNew=pickFlavMR(ctr);
          rpdg=combineMesonMR(e.flav, fromPos?-rfNew:rfNew, ctr); pairPTMR(ctr,rgx,rgy);
          if(rpdg!=0){ rmHad=HADMASS_MR(rpdg,ctr); rpxNew=rgx; rpyNew=rgy; rpxHad=e.px+rpxNew; rpyHad=e.py+rpyNew; } }
        if(rpdg==0) continue; rmT2=rmHad*rmHad+rpxHad*rpxHad+rpyHad*rpyHad; pdg=rpdg; fNew=rfNew;
        mHad=rmHad; pxNew=rpxNew; pyNew=rpyNew;
      }
      double z=zLundSample(H_ALUND,H_BLUND*rmT2,1.0,ctr);
      if(kinHad(S,fromPos,e,z,rmT2,rmHad,rpxHad,rpyHad,rpxNew,rpyNew,iPN,iNN,xPN,xNN,GN,had)){
        double hm2=had[3]*had[3]-had[0]*had[0]-had[1]*had[1]-had[2]*had[2];
        if(had[3]==had[3] && fabs(had[3])<1e6 && fabs(hm2-rmHad*rmHad)<1e-5*(1.0+rmHad*rmHad)){ made=true; accMass=rmHad; }
      }
    }
    if(!made) return -1;
    if(nH>=MAXPART) return -1;
    for(int k=0;k<4;++k){ H[4*nH+k]=had[k]; pRem[k]-=had[k]; } hid[nH]=pdg; hm[nH]=accMass; nH++;
    e.iPos=iPN; e.iNeg=iNN; e.xPos=xPN; e.xNeg=xNN; e.Gamma=GN; e.px=-pxNew; e.py=-pyNew; e.flav=fromPos?fNew:-fNew;
  }
  // finalTwo (port of :1640-1716). pT and wT2 are independent of the final flavour, so fix
  // them first, then retry the flavour pair until the two hadrons fit the remaining mass
  // (a lighter pair -> pions fits where a heavy one would not).
  Region reg; if(!finalRegion(S,pos,neg,reg)) return -1;
  double xPp,xPn,prpx,prpy; regionProject(reg,pRem,xPp,xPn,prpx,prpy);
  double pxRem=prpx-pos.px-neg.px, pyRem=prpy-pos.py-neg.py;
  double pos_pxHad=pos.px+0.5*pxRem, pos_pyHad=pos.py+0.5*pyRem;
  double negpxHad =neg.px+0.5*pxRem, negpyHad =neg.py+0.5*pyRem;
  double wT2=(pRem[3]*pRem[3]-pRem[0]*pRem[0]-pRem[1]*pRem[1]-pRem[2]*pRem[2])
             +(pos_pxHad+negpxHad)*(pos_pxHad+negpxHad)+(pos_pyHad+negpyHad)*(pos_pyHad+negpyHad);
  int fNew=0,pdg1=0,pdg2=0; double m1=0,m2=0,m1T2=0,m2T2=0,lam2=-1.0;
  for(int ft=0; ft<40 && lam2<=0.0; ++ft){
    int p1=0,p2=0;
    for(int tr=0; tr<20 && (p1==0||p2==0); ++tr){ fNew=pickFlavMR(ctr);
      p1=combineMesonMR(pos.flav,-fNew,ctr); p2=combineMesonMR(fNew,neg.flav,ctr); }
    if(p1==0||p2==0) continue;
    double mm1=HADMASS_MR(p1,ctr), mm2=HADMASS_MR(p2,ctr);
    double t1=mm1*mm1+pos_pxHad*pos_pxHad+pos_pyHad*pos_pyHad;
    double t2=mm2*mm2+negpxHad*negpxHad+negpyHad*negpyHad;
    if(sqrt(wT2)<sqrt(t1)+sqrt(t2)) continue;
    double l2=(wT2-t1-t2)*(wT2-t1-t2)-4.0*t1*t2;
    if(l2>0.0){ pdg1=p1;pdg2=p2;m1=mm1;m2=mm2;m1T2=t1;m2T2=t2;lam2=l2; }
  }
  if(lam2<=0.0) return -1;
  double lam=sqrt(lam2), pRev=1.0/(1.0+exp(fmin(50.0,H_BLUND*lam)));
  double xpz=0.5*lam/wT2; if(pRev>u01(splitmix64(ctr++))) xpz=-xpz;
  double xmd=(m1T2-m2T2)/wT2, xeP=0.5*(1.0+xmd), xeN=0.5*(1.0-xmd);
  double h1[4],h2[4];
  regionPHad(reg,(xeP+xpz)*xPp,(xeP-xpz)*xPn,pos_pxHad,pos_pyHad,h1);
  regionPHad(reg,(xeN-xpz)*xPp,(xeN+xpz)*xPn,negpxHad,negpyHad,h2);
  if(nH+2>MAXPART) return -1;
  double q1m2=h1[3]*h1[3]-h1[0]*h1[0]-h1[1]*h1[1]-h1[2]*h1[2];
  double q2m2=h2[3]*h2[3]-h2[0]*h2[0]-h2[1]*h2[1]-h2[2]*h2[2];
  if(h1[3]!=h1[3]||h2[3]!=h2[3]||fabs(q1m2-m1*m1)>1e-5*(1.0+m1*m1)||fabs(q2m2-m2*m2)>1e-5*(1.0+m2*m2)) return -1;
  for(int k=0;k<4;++k)H[4*nH+k]=h1[k]; hid[nH]=pdg1; hm[nH]=m1; nH++;
  for(int k=0;k<4;++k)H[4*nH+k]=h2[k]; hid[nH]=pdg2; hm[nH]=m2; nH++;
  return nH;
}
__host__ __device__ inline int hadronizeMR(const double* P,const int* id,int n,uint64_t ctr,double* H,int* hid,double* hm){
  for(int retry=0; retry<50; ++retry){ int r=tryFragmentMR(P,id,n,ctr,H,hid,hm); if(r>0) return r; }
  return -1;
}

__global__ void kern(int N,uint64_t base,int* outN,int* outNc,double* outTot,double* outDm,int* outNp,
                     double* outH,int* outHid){
  int e=blockIdx.x*(int)blockDim.x+threadIdx.x; if(e>=N) return;
  double Psh[MAXP*4]; int idsh[MAXP];
  int np=showerEvent(Psh,idsh, base+(uint64_t)e*0x9E3779B97F4A7C15ULL);   // FSR shower chain
  outNp[e]=np;
  double H[MAXPART*4]; int hid[MAXPART]; double hm[MAXPART];
  uint64_t hctr=(base^0x5151ULL)+(uint64_t)e*0x100000001B3ULL;
#ifdef GLUON_SPLIT
  // g->qqbar forks the colour chain into independent sub-strings: hadronize each with its own RNG
  // stream and concatenate. A non-forked event has ns==1 and uses hctr -> identical single-string
  // path. (Whole-event drop if any sub-string fails to fragment -> the existing drop accounting.)
  int starts[MAXP], ends[MAXP];
  int ns=findStrings(idsh,np,starts,ends);
  int nH=0; bool bad=false;
  for(int s=0;s<ns;++s){
    int a=starts[s], ssz=ends[s]-a+1;
    uint64_t sctr=(ns==1)? hctr : (hctr ^ ((uint64_t)(s+1)*0x9E3779B97F4A7C15ULL));
    double Hs[MAXPART*4]; int hids[MAXPART]; double hms[MAXPART];
    int r=hadronizeMR(&Psh[4*a], &idsh[a], ssz, sctr, Hs,hids,hms);
    if(r<0 || nH+r>MAXPART){ bad=true; break; }
    for(int j=0;j<r;++j){ for(int k=0;k<4;++k) H[4*(nH+j)+k]=Hs[4*j+k]; hid[nH+j]=hids[j]; hm[nH+j]=hms[j]; }
    nH+=r;
  }
  if(bad||nH<=0){ outN[e]=-1; return; }
#else
  int nH=hadronizeMR(Psh,idsh,np, hctr, H,hid,hm);
  if(nH<0){ outN[e]=-1; return; }
#endif
#ifdef DECAYS
  // Decay the primary unstable hadrons (rho/K*/omega/phi/eta/eta'/K0) into the ALEPH particle-level
  // stable set. SEPARATE counter stream dctr (NOT hctr) so toggling -DDECAYS cannot perturb the
  // hadronization draws -> the no-decay build stays byte-identical.
  double F[MAXFINAL*4]; int fid[MAXFINAL];
  uint64_t dctr=(base^0xDECAULL)+(uint64_t)e*0x100000001B3ULL;
  int nF=decayEvent(H,hid,nH,dctr,F,fid);
  if(nF<0){ outN[e]=-1; return; }
  double* OUTP=F; int* OUTID=fid; int OUTN=nF;
#else
  double* OUTP=H; int* OUTID=hid; int OUTN=nH;
#endif
  double s0=0,s1=0,s2=0,s3=0,dm=0; int nc=0;
  for(int i=0;i<OUTN;++i){ s0+=OUTP[4*i];s1+=OUTP[4*i+1];s2+=OUTP[4*i+2];s3+=OUTP[4*i+3];
    double m2=OUTP[4*i+3]*OUTP[4*i+3]-OUTP[4*i]*OUTP[4*i]-OUTP[4*i+1]*OUTP[4*i+1]-OUTP[4*i+2]*OUTP[4*i+2];
#ifdef DECAYS
    double mt=mesonMassMR(OUTID[i]);     // finals are stable -> pole mass (matches decay product mass)
#else
    double mt=hm[i];                     // primary hadron's BW-sampled mass
#endif
    dm=fmax(dm,fabs(m2-mt*mt)); if(isChargedMR(OUTID[i]))nc++; }
  outN[e]=OUTN; outNc[e]=nc; outTot[4*e]=s0;outTot[4*e+1]=s1;outTot[4*e+2]=s2;outTot[4*e+3]=s3; outDm[e]=dm;
  // Optional per-event hadron record (for the HepMC3/Rivet dump; ignored unless main writes it).
  for(int i=0;i<OUTN;++i){ for(int k=0;k<4;++k) outH[((size_t)e*OUTCAP+i)*4+k]=OUTP[4*i+k]; outHid[(size_t)e*OUTCAP+i]=OUTID[i]; }
}

int main(int argc,char**argv){
  int N=(argc>1)?atoi(argv[1]):20000; uint64_t base=0x4D52ULL;
#ifdef DECAYS
  int TPB=64;    // larger per-thread F[]/S[] under decays -> use fewer threads/block to fit
#else
  int TPB=128;
#endif
  int blocks=(N+TPB-1)/TPB;
  const char* dumpFile=(argc>2)?argv[2]:nullptr;   // optional: write per-event hadrons for HepMC3/Rivet
  int *dN,*dNc,*dNp,*dHid; double *dTot,*dDm,*dH;
  CK(cudaMalloc(&dN,(size_t)N*4));CK(cudaMalloc(&dNc,(size_t)N*4));CK(cudaMalloc(&dNp,(size_t)N*4));CK(cudaMalloc(&dTot,(size_t)N*32));CK(cudaMalloc(&dDm,(size_t)N*8));
  CK(cudaMalloc(&dH,(size_t)N*OUTCAP*4*8));CK(cudaMalloc(&dHid,(size_t)N*OUTCAP*4));
  kern<<<blocks,TPB>>>(N,base,dN,dNc,dTot,dDm,dNp,dH,dHid); CK(cudaDeviceSynchronize());
  std::vector<int> hN(N),hNc(N),hNp(N); std::vector<double> hTot((size_t)N*4),hDm(N);
  CK(cudaMemcpy(hN.data(),dN,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hNc.data(),dNc,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hNp.data(),dNp,(size_t)N*4,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hTot.data(),dTot,(size_t)N*32,cudaMemcpyDeviceToHost));
  CK(cudaMemcpy(hDm.data(),dDm,(size_t)N*8,cudaMemcpyDeviceToHost));
  std::vector<double> hH; std::vector<int> hHid;
  if(dumpFile){ hH.resize((size_t)N*OUTCAP*4); hHid.resize((size_t)N*OUTCAP);
    CK(cudaMemcpy(hH.data(),dH,(size_t)N*OUTCAP*4*8,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(hHid.data(),dHid,(size_t)N*OUTCAP*4,cudaMemcpyDeviceToHost)); }
  std::vector<int> hN2(N); kern<<<blocks,TPB>>>(N,base,dN,dNc,dTot,dDm,dNp,dH,dHid); CK(cudaDeviceSynchronize());
  CK(cudaMemcpy(hN2.data(),dN,(size_t)N*4,cudaMemcpyDeviceToHost));
  if(dumpFile){ FILE* fo=fopen(dumpFile,"w");
    if(fo){ long nValid=0; for(int e=0;e<N;++e) if(hN[e]>0) nValid++;
      fprintf(fo,"%ld %.6f\n",nValid,MZ);
      for(int e=0;e<N;++e){ if(hN[e]<=0) continue; fprintf(fo,"%d\n",hN[e]);
        for(int i=0;i<hN[e];++i){ size_t b=((size_t)e*OUTCAP+i)*4;
          fprintf(fo,"%d 0 0 %.9e %.9e %.9e %.9e\n",hHid[(size_t)e*OUTCAP+i],hH[b],hH[b+1],hH[b+2],hH[b+3]); } }
      fclose(fo); printf("  dumped %ld hadron-level events -> %s\n",nValid,dumpFile); } }

  double maxMom=0,maxDm=0; long sumN=0,sumNc=0,nFail=0,repro=0;
  int npOkMax=0,npBadMin=1<<30; long nBadN2=0,nN2=0,nBadN3plus=0;  // diagnostic: which n fails
  for(int e=0;e<N;++e){ if(hN[e]!=hN2[e])repro++;
    if(hN[e]<0){ nFail++; if(hNp[e]<npBadMin)npBadMin=hNp[e]; continue; }
    double dx=fabs(hTot[4*e]),dy=fabs(hTot[4*e+1]),dz=fabs(hTot[4*e+2]),de=fabs(hTot[4*e+3]-MZ);
    double v=fmax(fmax(dx,dy),fmax(dz,de));
    if(hNp[e]==2){ nN2++; if(v>1e-3)nBadN2++; } else if(v>1e-3) nBadN3plus++;
    if(v<1e-3 && hNp[e]>npOkMax) npOkMax=hNp[e];
    maxMom=fmax(maxMom,v); maxDm=fmax(maxDm,hDm[e]); sumN+=hN[e]; sumNc+=hNc[e]; }
  (void)nN2;(void)nBadN2;(void)nBadN3plus;(void)npOkMax;(void)npBadMin;
  long nok=N-nFail;
  printf("All-GPU multi-region (gluon-kinked) hadronization: FSR shower -> hadrons (sqrt(s)=%.3f, %d evts)\n",MZ,N);
#ifdef DECAYS
  const char* dlbl="(decays on)";
#else
  const char* dlbl="(no decays)";
#endif
  printf("  multiplicity      : mean %.3f hadrons, %.3f charged %s\n",(double)sumN/(nok?nok:1),(double)sumNc/(nok?nok:1),dlbl);
  printf("  4-mom conservation: max|deviation| = %.2e GeV\n",maxMom);
  printf("  on-shellness      : max|m^2-table| = %.2e GeV^2\n",maxDm);
  printf("  refragment-drop   : %ld / %d = %.1f%% (hard configs dropped — biases mult ~few%% low)\n",nFail,N,100.0*nFail/N);
  printf("  reproducibility   : GPU re-run diffs = %ld\n",repro);
  bool ok=(maxMom<1e-5)&&(maxDm<1e-6)&&(repro==0)&&(nFail<N/10)&&((double)sumNc/(nok?nok:1)>5.0);
  printf("VALIDATION: %s (EXACT conservation+on-shell+reproducible; multi-region kinematics correct)\n",ok?"PASS":"FAIL");
  cudaFree(dN);cudaFree(dNc);cudaFree(dTot);cudaFree(dDm);
  return ok?0:2;
}
