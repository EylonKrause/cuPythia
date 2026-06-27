// Pythia 8.317 reference for the cuPythia GPU FSR shower: e+e- -> Z -> q qbar at the
// Z pole, FINAL-STATE radiation ONLY (ISR/MPI/hadronization off), matched as closely as
// possible to cuPythia/pipeline/shower_fsr.cu so the thrust distributions can be compared.
//
// Matching choices: pTmin = 0.5 GeV, alphaSvalue = 0.1365, alphaSorder = 1 (Pythia
// defaults), Z -> u,d,s only (near-massless, like the GPU shower's massless partons),
// QED/weak FSR off. Honest, documented differences vs the GPU shower (shower_fsr.cu):
//  - Pythia keeps g->qqbar ON (nGluonToQuark default); the GPU shower omits it. We CANNOT
//    set nGluonToQuark=0 here: it sends SimpleTimeShower into an infinite loop in 8.317.
//  - Pythia matches alpha_s across flavour thresholds; the GPU shower fixes n_f=5.
//  - Pythia applies ME corrections to the first emission; the GPU shower does not.
// These are the residuals the thrust comparison exposes (and what closing the "fit" means).
//
// Build: g++ thrust_pythia.cc -o thrust_pythia $(../../pythia8317/bin/pythia8-config --cxxflags --libs)
// Run:   ./thrust_pythia [nEvents=200000]

#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;

int main(int argc,char**argv){
  int nEvt=(argc>1)?atoi(argv[1]):200000;
  bool mecOff=false, cmw=false;
  for(int a=2;a<argc;++a){ std::string s(argv[a]); if(s=="mecoff")mecOff=true; if(s=="cmw")cmw=true; }
  Pythia pythia;
  pythia.readString("Beams:idA = -11");                  // e+   (matches Pythia example main111)
  pythia.readString("Beams:idB = 11");                   // e-
  pythia.readString("Beams:eCM = 91.1876");
  pythia.readString("PDF:lepton = off");                 // e+- pointlike, no QED ISR
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");   // e+e- -> gamma*/Z
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3");               // Z -> u,d,s (near-massless)
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = off");            // stop at the parton level
  pythia.readString("TimeShower:QEDshowerByQ = off");    // QCD FSR only
  pythia.readString("TimeShower:pTmin = 0.5");            // match the GPU shower cutoff
  if(mecOff) pythia.readString("TimeShower:MEcorrections = off"); // LL shower, like the GPU shower
  if(cmw)    pythia.readString("TimeShower:alphaSuseCMW = on");   // soft-coherence NLL rescaling
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");
  if(!pythia.init()){ printf("Pythia init failed\n"); return 1; }

  const int NB=20; const double TMAX=0.5; long hist[NB]={0}; long nAcc=0; double sum1mT=0;
  Thrust thr(2);                                          // select=2: all final particles
  for(int e=0;e<nEvt;++e){ if(!pythia.next()) continue;
    if(!thr.analyze(pythia.event)) continue;
    double omt=1.0-thr.thrust(); sum1mT+=omt; nAcc++;
    int b=(int)(omt/TMAX*NB); if(b<0)b=0; if(b>=NB)b=NB-1; hist[b]++; }

  const char* fn = cmw ? "thrust_pythia_mecoff_cmw.dat"
                       : (mecOff ? "thrust_pythia_mecoff.dat" : "thrust_pythia.dat");
  FILE* f=fopen(fn,"w");
  if(f){ fprintf(f,"# (1-T)_low  (1-T)_high  normalised_density   [Pythia 8.317 FSR, MEcorr=%s, %ld evts]\n",
            mecOff?"off":"on", nAcc);
    for(int b=0;b<NB;++b) fprintf(f,"%.4f %.4f %.6e\n",b*TMAX/NB,(b+1)*TMAX/NB,hist[b]/((double)nAcc*(TMAX/NB)));
    fclose(f); }
  printf("Pythia 8.317 FSR reference (MEcorr=%s): %ld events, <1-T> = %.4f  -> %s\n",
         mecOff?"off":"on", nAcc, sum1mT/nAcc, fn);
  return 0;
}
