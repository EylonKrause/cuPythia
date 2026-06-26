// main303.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Peter Skands <peter.skands@monash.edu>

// Keywords: electron-positron; histograms; openmp

// Simple example of final-state showers at LEP 1, comparing various options
// for higher-order corrections (non-singular 1st-order and 2nd-order terms
// in the branching kernels, with and without 1st-order MECs, and different
// options for IR regularisations of alphaS) to a few reference settings.
//
// Exploits option for parallel runs using OpenMP if the --with-openmp option
// was used when running the Pythia configure script.

// 1st-order non-singular terms.
const double cEmitQ = 2.0;
const double cEmitG = 2.0;
const double cSplit = 0.0;
// 2nd-order splitting-function terms.
const double hEmitHard  = 0.0;
const double hEmitColl  = 0.0;
const double hEmitSoft  = 0.0;
const double hSplitHard = 0.0;
const double hSplitColl = 0.0;
// ME corrections on/off.
const bool   MEcorrections = true;
// IR cutoff of shower.
const double cutoff = 1.0;
// alphaS(mZ) value, running order, and scheme choice.
const double alphaSmZ     = 0.12;
const int    alphaSorder  = 2;
const bool   alphaSuseCMW = true;
// Maximum alphaS value (in the IR)
const double alphaSmax = 4.0;
// Renormalization-scale shift factor (in units of Lambda3).
const double alphaSrenormShift = 1.0;

// Include Pythia8 header(s) and namespace.
#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Main program.
int main() {

  // Arrays to store averages, outside of runs.
  const int NRUNS = 3;
  double avgNch[NRUNS], avgTau[NRUNS], avgMaj[NRUNS], avgMin[NRUNS],
    avgObl[NRUNS];

  // Loop over generator runs. Use OpenMP parallelisation if enabled.
#ifdef OPENMP
#pragma omp parallel for
#endif
  for (int iRun = 0; iRun < NRUNS; ++iRun) {

    // Define Pythia 8 generator
    Pythia pythia;

    // Settings for hadronic Z decays.
    pythia.readString("Beams:idA =  11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
    pythia.readString("PDF:lepton = off");
    pythia.readString("SpaceShower:QEDshowerByL = off");
    pythia.readString("Beams:eCM = 91.2");

    //----------------------------------------------------------------------

    // Shorthands.
    Event& event = pythia.event;
    Settings& settings = pythia.settings;

    // Extract settings to be used in the main program.
    int    nEvent     = 100000;
    int    nAbort     = settings.mode("Main:timesAllowErrors");

    if ( iRun == 0 ) {
      // Run 0: Default (Monash) settings.
    } else if ( iRun == 1 ) {
      // Run 1: alphaS = 0.118, 2-loop running, and CMW.
      pythia.readString("TimeShower:alphaSvalue = 0.118");
      pythia.readString("TimeShower:alphaSorder = 2");
      pythia.readString("TimeShower:alphaSuseCMW = on");
    } else {
      // User-specified settings for effective splitting kernels.
      settings.parm("TimeShower:cEmitQ", cEmitQ);
      settings.parm("TimeShower:cEmitG", cEmitG);
      settings.parm("TimeShower:cSplit", cSplit);
      settings.parm("TimeShower:hEmitHard", hEmitHard);
      settings.parm("TimeShower:hEmitColl", hEmitColl);
      settings.parm("TimeShower:hEmitSoft", hEmitSoft);
      settings.parm("TimeShower:hSplitHard", hSplitHard);
      settings.parm("TimeShower:hSplitColl", hSplitColl);
      settings.flag("TimeShower:MEcorrections", MEcorrections);
      // Shower IR cutoff.
      settings.parm("TimeShower:pTmin", cutoff);
      // User-specified settings for alphaS.
      settings.parm("TimeShower:alphaSvalue", alphaSmZ);
      settings.mode("TimeShower:alphaSorder", alphaSorder);
      settings.flag("TimeShower:alphaSuseCMW", alphaSuseCMW);
      settings.parm("TimeShower:alphaSmax", alphaSmax);
      settings.parm("TimeShower:alphaSrenormShift", alphaSrenormShift);
    }

    // Prevent multithreadhing during init (for nicer output).
    bool ok = true;
#ifdef OPENMP
#pragma omp critical
#endif
    {
      // If Pythia fails to initialize, skip this run.
      cout << "Initialising PYTHIA for Run "<< iRun << endl;
      ok = pythia.init();
      cout << "Finished PYTHIA Initialisation for Run "<< iRun << endl;
    }
    if (!ok) continue;

    //----------------------------------------------------------------------

    // Define instance of Thrust for event analysis.
    Thrust thrust(1);

    //----------------------------------------------------------------------

    // Define and fill a histogram for effective alphaS value.

    // AlphaS settings (for plotting).
    AlphaStrong alpS;
    alpS.init( alphaSmZ, alphaSorder, 6, alphaSuseCMW,
      alphaSmax, alphaSrenormShift );

    // Plot value of effective AlphaS, logarithmic in Q.
    int nBinsAS = 100;
    double qMin = 0.1;
    double qMax = 100.;
    Hist hAlphaS("alphaS", nBinsAS, qMin, qMax, true);
    for (int iBin = 1; iBin <= nBinsAS; ++iBin) {
      double q  = hAlphaS.getBinCenter( iBin );
      double aS = alpS.alphaS( pow2(q) );
      hAlphaS.fill( q, aS );
    }

    // Don't mix threads when printing the following.
#ifdef OPENMP
#pragma omp critical
#endif
    {
      cout << hAlphaS;
    }

    //----------------------------------------------------------------------

    // EVENT GENERATION LOOP.
    // Generation, event-by-event printout, analysis, and histogramming.

    // Counters.
    double weight = 1.0;
    double sumWeights = 0.0;
    double sumTau(0), sumMaj(0), sumMin(0), sumObl(0), sumNch(0);

    // Begin event loop.
    int iAbort = 0;
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      bool aborted = !pythia.next();
      if (aborted){
        event.list();
        if (++iAbort < nAbort) {
          continue;
        }
        cout << " Event generation aborted prematurely, owing to error!\n";
        cout << " Event number was : " << iEvent << endl;
        break;
      }

      // Check for weights.
      weight = pythia.info.weight();
      sumWeights += weight;

      // Count FS charged hadrons, partons, and quarks.
      int nCharged  = 0;
      for (int i = 1; i<event.size(); i++) {
        // Final-state partons.
        if ( !event[i].isFinal() ) continue;
        if ( event[i].isCharged() ) nCharged++;
      }

      // Fill charged multiplicity for this run.
      sumNch += nCharged * weight;

      // Thrust analysis.
      thrust.analyze( event );
      sumTau += weight * (1. - thrust.thrust());
      sumMaj += weight * thrust.tMajor();
      sumMin += weight * thrust.tMinor();
      sumObl += weight * thrust.oblateness();

    }

    //----------------------------------------------------------------------

    // POST-RUN FINALIZATION.
    // Normalization, Statistics, Output.

    // Normalize histograms to effective number of positive-weight events.
    double normFac = 1.0/sumWeights;

    // Don't mix threads when saving to arrays.
#ifdef OPENMP
#pragma omp critical
#endif
    {
      avgNch[ iRun ] = sumNch * normFac;
      avgTau[ iRun ] = sumTau * normFac;
      avgMaj[ iRun ] = sumMaj * normFac;
      avgMin[ iRun ] = sumMin * normFac;
      avgObl[ iRun ] = sumObl * normFac;
    }
  }

  // Wait for all threads to complete before printing summary:
#ifdef OPENMP
#pragma omp barrier
#endif

  // Write out averages and RMSD widths.
  const int nChar = 10;
  cout << "\n Summary of runs: " << endl;
  for (int iRun = 0; iRun < NRUNS; ++iRun) {
    cout << "  iRun = " << num2str(iRun,2) <<": "
         << "  <tau> = " << num2str(avgTau[iRun], nChar)
         << "  <maj> = " << num2str(avgMaj[iRun], nChar)
         << "  <min> = " << num2str(avgMin[iRun], nChar)
         << "  <obl> = " << num2str(avgObl[iRun], nChar)
         << "  <nCh> = " << num2str(avgNch[iRun], nChar)
         << endl;
  }

  cout << endl;
  // Done.
  return 0;
}
