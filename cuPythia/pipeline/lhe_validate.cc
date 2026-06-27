// Validate the cuPythia LHE output (lhe_writer.cu) by reading it straight back into
// Pythia 8.317 (Beams:frameType=4) and showering+hadronizing every event. This proves the
// file is spec-valid LHEF AND physically usable: Pythia must parse <init>/<event>, accept
// the colour flow, and conserve total 4-momentum (= 2*Ebeam) for every event.
//
// Build: g++ lhe_validate.cc -o lhe_validate $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
// Run:   ./lhe_validate [events.lhe]

#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cmath>
using namespace Pythia8;

int main(int argc,char**argv){
  const char* path=(argc>1)?argv[1]:"events.lhe";
  Pythia pythia;
  pythia.readString("Beams:frameType = 4");
  pythia.settings.word("Beams:LHEF", path);
  pythia.readString("PartonLevel:MPI = off");   // keep it light; not needed for the I/O check
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");
  if(!pythia.init()){ printf("LHE VALIDATE: FAIL (Pythia could not init from %s)\n",path); return 2; }

  double Etot=2.0*6500.0;       // matches the LHE <init> beam energies
  long nOk=0, nAbort=0; double maxImb=0;
  for(int e=0; ; ++e){
    if(!pythia.next()){
      if(pythia.info.atEndOfFile()) break;
      if(++nAbort > 10) break;   // tolerate a few, but a valid file shouldn't abort
      continue;
    }
    nOk++;
    double s0=0,s1=0,s2=0,s3=0;
    for(int i=0;i<pythia.event.size();++i) if(pythia.event[i].isFinal()){
      s0+=pythia.event[i].px(); s1+=pythia.event[i].py();
      s2+=pythia.event[i].pz(); s3+=pythia.event[i].e(); }
    double imb=fabs(s0)+fabs(s1)+fabs(s2)+fabs(s3-Etot);
    if(imb>maxImb) maxImb=imb;
  }
  printf("LHE readback via Pythia 8.317 (frameType=4): %s\n",path);
  printf("  events read+showered+hadronized = %ld   aborted = %ld\n", nOk, nAbort);
  printf("  max total 4-momentum imbalance  = %.2e GeV (expect ~ MC/remnant level)\n", maxImb);
  bool ok = (nOk>0) && (nAbort==0) && (maxImb < 1.0);
  printf("LHE VALIDATE: %s (spec-valid + colour-valid + showerable + conserving)\n", ok?"PASS":"FAIL");
  return ok?0:2;
}
