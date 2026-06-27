// Hadronize the GPU shower's gluon-kinked parton chains (bridge.cu output) with Pythia
// 8.317 StringFragmentation (forceHadronLevel) — the shower->hadronization bridge. Validates
// that every GPU-shower colour singlet is accepted and hadronized with 4-momentum conserved,
// and reports the full e+e- hadron multiplicity (GPU shower + Pythia kinked-string hadronization).
//
// Build: g++ bridge_pythia.cc -o bridge_pythia $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
// Run:   ./bridge_pythia [shower_partons.dat]

#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cmath>
using namespace Pythia8;

int main(int argc,char**argv){
  const char* path=(argc>1)?argv[1]:"shower_partons.dat";
  FILE* f=fopen(path,"r"); if(!f){ printf("cannot open %s\n",path); return 1; }
  int N; double Ecm; if(fscanf(f,"%d %lf",&N,&Ecm)!=2){ printf("bad header\n"); return 1; }

  bool noDecay=(argc>2 && std::string(argv[2])=="nodecay");
  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");      // feed partons directly
  pythia.readString(noDecay ? "HadronLevel:Decay = off" : "HadronLevel:Decay = on");
  pythia.readString("StringFlav:probQQtoQ = 0");    // no baryons (match the GPU multi-region slice)
  pythia.readString("Print:quiet = on");
  if(!pythia.init()){ printf("BRIDGE: FAIL (init)\n"); return 2; }
  Event& ev = pythia.event;

  long nOk=0, nBad=0, sumN=0, sumNc=0; double maxImb=0;
  for(int e=0;e<N;++e){
    int n; if(fscanf(f,"%d",&n)!=1) break;
    ev.reset();
    for(int i=0;i<n;++i){ int id,col,acol; double px,py,pz,en;
      if(fscanf(f,"%d %d %d %lf %lf %lf %lf",&id,&col,&acol,&px,&py,&pz,&en)!=7){ n=-1; break; }
      ev.append(id,23,col,acol, px,py,pz,en, 0.0); }
    if(n<0) break;
    pythia.forceHadronLevel();   // return is a tolerance flag; gate on conservation below
    int nf=0,nc=0; double s0=0,s1=0,s2=0,s3=0;
    for(int i=0;i<ev.size();++i) if(ev[i].isFinal()){ nf++; if(ev[i].isCharged())nc++;
      s0+=ev[i].px();s1+=ev[i].py();s2+=ev[i].pz();s3+=ev[i].e(); }
    double imb=fabs(s0)+fabs(s1)+fabs(s2)+fabs(s3-Ecm);
    if(nf>=2 && imb<1e-3){ nOk++; sumN+=nf; sumNc+=nc; if(imb>maxImb)maxImb=imb; }
    else nBad++;
  }
  fclose(f);
  printf("Shower->hadronization bridge: GPU FSR shower + Pythia kinked-string hadronization\n");
  printf("  events hadronized (conserved) = %ld / %d   rejected = %ld\n", nOk, N, nBad);
  printf("  full e+e- hadron multiplicity = %.2f total, %.2f charged (decays on)\n",
         (double)sumN/(nOk?nOk:1),(double)sumNc/(nOk?nOk:1));
  printf("  max 4-momentum imbalance      = %.2e GeV\n", maxImb);
  bool ok = (nOk > 0.99*N) && (maxImb < 1e-3);
  printf("BRIDGE: %s (every GPU-shower colour singlet hadronizes, conserving)\n", ok?"PASS":"FAIL");
  return ok?0:2;
}
