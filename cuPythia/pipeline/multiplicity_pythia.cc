// Pythia 8.317 reference for the cuPythia GPU hadronization slice (hadronize.cu).
// Hadronizes the SAME object — a single straight u-ubar colour singlet at fixed sqrt(s),
// hadron level only — and reports primary-hadron multiplicity for an apples-to-apples
// (modulo listed residuals) comparison.
//
// Matched config: ProcessLevel/PartonLevel OFF (we feed the partons), HadronLevel:Decay OFF
// (primary hadrons only, like the GPU slice), StringFlav:probQQtoQ=0 (no baryons, like the
// slice), default fragmentation params (pinned). Documented residuals vs the slice: Pythia
// gives vector mesons Breit-Wigner masses (the slice uses pole masses); Pythia's full
// StringFlav includes the (default-zero-rate) L=1 multiplets the slice omits.
//
// Build: g++ multiplicity_pythia.cc -o multiplicity_pythia $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
// Run:   ./multiplicity_pythia [nEvents=50000] [sqrtS=91.1876]

#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;

int main(int argc,char**argv){
  int nEvt=(argc>1)?atoi(argv[1]):50000;
  double rootS=(argc>2)?atof(argv[2]):91.1876, E=0.5*rootS;
  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");      // feed partons directly
  pythia.readString("HadronLevel:Decay = off");     // primary hadrons only
  pythia.readString("StringFlav:probQQtoQ = 0");    // no diquarks -> no baryons
  pythia.readString("StringZ:aLund = 0.68");
  pythia.readString("StringZ:bLund = 0.98");
  pythia.readString("StringPT:sigma = 0.335");
  pythia.readString("StringFlav:probStoUD = 0.217");
  pythia.readString("StringFlav:mesonUDvector = 0.50");
  pythia.readString("StringFlav:mesonSvector = 0.55");
  pythia.readString("Print:quiet = on");
  if(!pythia.init()){ printf("init failed\n"); return 1; }

  Event& event = pythia.event;
  long sumN=0,sumNc=0,nAcc=0,nSkip=0;
  for(int e=0;e<nEvt;++e){
    event.reset();
    // u (col 101) along +z, ubar (acol 101) along -z — a colour singlet, massless ends.
    event.append( 2,23,101,0,  0.0,0.0, E,E, 0.0);
    event.append(-2,23,0,101,  0.0,0.0,-E,E, 0.0);
    pythia.forceHadronLevel();   // return is an internal tolerance flag; gate on conservation below
    int n=0,nc=0; double s0=0,s1=0,s2=0,s3=0;
    for(int i=0;i<event.size();++i) if(event[i].isFinal()){
      n++; if(event[i].isCharged()) nc++;
      s0+=event[i].px(); s1+=event[i].py(); s2+=event[i].pz(); s3+=event[i].e(); }
    double dev=fabs(s0)+fabs(s1)+fabs(s2)+fabs(s3-rootS);
    if(n<2 || dev>0.1){ nSkip++; continue; }    // only count cleanly-conserved hadronizations
    sumN+=n; sumNc+=nc; nAcc++;
  }
  printf("  (%ld events skipped on conservation > 0.1 GeV)\n",nSkip);
  printf("Pythia 8.317 single u-ubar string hadronization (sqrt(s)=%.4f GeV, %ld events)\n",rootS,nAcc);
  printf("  primary multiplicity: mean %.3f hadrons, %.3f charged\n",
         (double)sumN/nAcc,(double)sumNc/nAcc);
  return 0;
}
