// main372.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjörn Sjöstrand <torbjorn.sjostrand@fysik.lu.se>
//          Christian T. Preuss <christian.preuss@uni-goettingen.de>

// Keywords: top; toponium; spin correlations

// Study spin correlations of the t-tbar system in t-tbar production
// below threshold.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

int main() {

  // Number of events. Full event simulation or simplified.
  int  nEvent = 1000000;
  bool fullEvents = false;

  // Book histograms.
  Hist coslmlpOff[6], coslmlpOn[6];
  for (int eBin = 0; eBin < 6; ++eBin) {
    coslmlpOff[eBin].book( "cosTheta(l-, l+), off", 40, -1., 1.);
    coslmlpOn[eBin].book( "cosTheta(l-, l+), on", 40, -1., 1.);
  }

  // Loop over angular correlations on or off.
  for (int iOO = 0; iOO < 2; ++iOO) {

    // Off/on pseudoscalar angular correlations, and range where dampened out.
    bool   doPseudoscalar  = (iOO == 0) ? false : true;
    double eBegDamp        = 10.;
    double eEndDamp        = 20.;

    // Example of other possible setup values.
    int    topModel        = 2;
    double mt              = 172.5;
    double gammat          = 1.34;
    double gammatGreen     = 1.34;
    double thresholdRegion = 10.;
    int    alphasOrder     = 2;
    double alphasValue     = 0.118;
    double ggSingletFrac   = 2./7.;
    double qqSingletFrac   = 0.;

    // Generator.
    Pythia pythia;

    // Process setup: select sequence t -> b W -> b lepton antilepton.
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("6:onMode = off");
    pythia.readString("6:onIfAll = 24 5");
    pythia.readString("24:onMode = off");
    pythia.readString("24:onIfAny = 11 13");

    // Kinematics setup: threshold region only.
    pythia.readString("Beams:eCM = 13000.");
    pythia.readString("PhaseSpace:mHatMin = 300.");
    pythia.readString("PhaseSpace:mHatMax = 400.");

    // Top threshold model setup.
    pythia.settings.mode("TopThreshold:model", topModel);
    pythia.particleData.m0(6, mt);
    pythia.readString("6:doForceWidth = true");
    pythia.particleData.mWidth( 6, gammat);
    pythia.settings.parm("TopThreshold:tWidthGreen", gammatGreen);
    pythia.settings.parm("TopThreshold:thrRegion", thresholdRegion);
    pythia.settings.mode("TopThreshold:alphasOrder", alphasOrder);
    pythia.settings.parm("TopThreshold:alphasValue", alphasValue);
    pythia.settings.parm("TopThreshold:ggSingletFrac", ggSingletFrac);
    pythia.settings.parm("TopThreshold:qqSingletFrac", qqSingletFrac);
    pythia.settings.flag("TopThreshold:pseudoscalar", doPseudoscalar);
    pythia.settings.parm("TopThreshold:psEBeginDamp", eBegDamp);
    pythia.settings.parm("TopThreshold:psEEndDamp", eEndDamp);

    // Switch off (most) code parts not relevant here. Reduce printout.
    if (!fullEvents) {
      pythia.readString("PartonLevel:ISR = off");
      pythia.readString("PartonLevel:FSR = off");
      pythia.readString("PartonLevel:MPI = off");
      pythia.readString("HadronLevel:all = off");
    }
    pythia.readString("Next:numberCount = 100000");

    // If Pythia fails to initialize, exit with error.
    if (!pythia.init()) return 1;

    // Begin event loop. Generate event. Skip if error. List first one.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) continue;

      // Create working copy of process record and boost it longitudinally
      // to the t + tbar rest frame. List first boosted event.
      Event work = pythia.process;
      double betaz = (work[5].pz() + work[6].pz()) /
        (work[5].e() + work[6].e());
      work.bst(0., 0., -betaz);
      if (iEvent == 0) {
        cout << "\n Listing of work, i.e. process boosted to "
          "the t + tbar rest frame:";
        work.list();
      }

      // Order ttbar system and the decay products (treg = t).
      int itreg = (work[5].id() == 6) ? 5 : 6;
      int itbar = 11 - itreg;
      int iWpos = work[itreg].daughter1();
      int iWneg = work[itbar].daughter1();
      // The W+- decay products, ferm = fermion, anti = antifermion.
      int iWposferm = work[iWpos].daughter1();
      int iWposanti = work[iWpos].daughter2();
      if (work[iWposferm].id() < 0) swap( iWposferm, iWposanti);
      int iWnegferm = work[iWneg].daughter1();
      int iWneganti = work[iWneg].daughter2();
      if (work[iWnegferm].id() < 0) swap( iWnegferm, iWneganti);

      // Angle between l- and l+  in subprocess rest frame.
      double cosThetalmlp  = costheta( work[iWnegferm].p(),
        work[iWposanti].p() );

      // Energy relative to nominal threshold, and derived energy bins.
      double eThr = pythia.info.toponiumE;
      int eBin = 0;
      if (eThr > eEndDamp) eBin = 4;
      else if (eThr > 0.5 * (eBegDamp + eEndDamp)) eBin = 3;
      else if (eThr > eBegDamp) eBin = 2;
      else if (eThr > 0.) eBin = 1;

      // Fill histograms, with event weight just in case.
      double wt   = pythia.info.weight();
      if (iOO == 0) coslmlpOff[eBin].fill( cosThetalmlp, wt);
      else          coslmlpOn[eBin].fill( cosThetalmlp, wt);

    // End of event loop. Statistics. End of off/on pseudoscalar angles.
    }
    pythia.stat();
  }

  // Combine threshold region. Normalize histograms to unit area.
  coslmlpOff[5] = coslmlpOff[0] + coslmlpOff[1];
  coslmlpOn[5] = coslmlpOn[0] + coslmlpOn[1];
  for (int eBin = 0; eBin < 6; ++eBin) {
    coslmlpOff[eBin].normalize(20.);
    coslmlpOn[eBin].normalize(20.);
  }

  // Plot histograms.
  HistPlot hpl("plot372");
  string linest[5] = { "h,black", "h,red", "h,blue", "h,magenta", "h,olive"};
  string captst[5] = { "below threshold", "below damping region",
    "begin damping region", "end damping region", "above damping region"};

  // All histograms with off or on, respectively.
  hpl.frame("fig372",
    "cosine of opening angle between two leptons, no correlations",
    "$\\cos\\theta_{\\ell^+\\ell^-}$",  "normalized rate");
  for (int eBeg = 0; eBeg < 5; ++eBeg)
    hpl.add(coslmlpOff[eBeg], linest[eBeg], captst[eBeg]);
  hpl.plot( -1., 1., 0., 1.);
  hpl.frame("",
    "cosine of opening angle between two leptons, with correlations",
    "$\\cos\\theta_{\\ell^+\\ell^-}$",  "normalized rate");
  for (int eBeg = 0; eBeg < 5; ++eBeg)
    hpl.add(coslmlpOn[eBeg], linest[eBeg], captst[eBeg]);
  hpl.plot( -1., 1., 0., 1.);

  // Threshold region, off vs on.
  hpl.frame("",
    "cosine of opening angle between two leptons, without/with correlations",
    "$\\cos\\theta_{\\ell^+\\ell^-}$",  "normalized rate");
  hpl.add(coslmlpOff[5], "h,red", "threshold region, no correlations");
  hpl.add(coslmlpOn[5], "h,blue", "threshold region, with correlations");
  hpl.plot( -1, 1., 0., 1.);

  // Done.
  return 0;

}
