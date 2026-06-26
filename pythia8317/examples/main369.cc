// main369.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: top; toponium; openmp;

// Plot ttbar system properties contrasting five different scenarios:
// 0 : pure Born baseline.
// 1 : Coulomb enhancement.
// 2 : Green's function, full top width in Breit-Wigner and matrix element.
// 3 :  Green's function, reduced top width in Breit-Wigner.
// 4 :  Green's function, reduced top width in matrix elements.
// Plots all scenarios in one frame, which may be too much.

// Runs in parallel if Pythia configured with --with-openmp.
// The same analysis is also provided with main370 but using PythiaParallel.

// Based on threshold factors according to
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

// OpenMP includes an user-specifiable number of threads requested.
// (Will be capped by system maximum. Use -1 to use all available.)
#ifdef OPENMP
#include <omp.h>
const int nThreadsRequested = -1;
#endif

//============================================================================

int main() {

  // Restrict to threshold region only or not.
  bool thresholdOnly   = true;

  // Full event generation or only cross-section-oriented studies.
  bool fullEvents      = false;

  // Number of events.
  int nEvent           = 1000000;

  // LHC collision energy.
  double eCM           = 13000.;

  // Example of common possible setup values.
  double mt            = 172.5;
  double thrRegion     = 10.;
  int    alphasOrder   = 2;
  double alphasValue   = 0.118;
  double ggSingletFrac = 2. / 7.;
  double qqSingletFrac = 0.;

  // Histograms.
  Hist   mHatLow[5], mHatLowZoom[5], mThrOrig[5], mThrOrigZoom[5];
  for (int scenario = 0; scenario < 5; ++scenario) {
    mHatLow[scenario].book("ttbar invariant mass, low", 100, 300., 400.);
    mHatLowZoom[scenario].book("ttbar invariant mass, low", 100, 330., 355.);
    mThrOrig[scenario].book("original threshold", 100, -20., 80.);
    mThrOrigZoom[scenario].book("original threshold, zoom", 100, -10., 30.);
  }

  // Compare the five scenarios.
  // Automatically run in parallel if Pythia configured with flag
  // --with-openmp.
#ifdef OPENMP
  int nThreadsSys = omp_get_max_threads();
  int nThreadsMax = min(5, nThreadsSys);
  int nThreads    = nThreadsRequested;
  if (nThreads <= 0 || nThreads > nThreadsMax) nThreads = nThreadsMax;
  cout << " OMP parralelisation using "<<nThreads << " threads."
       << " System maximum = "<<nThreadsSys<<"." << endl;
#pragma omp parallel for num_threads(nThreads)
#endif
  for (int scenario = 0; scenario < 5; ++scenario) {

    // Scenario-specific settings choices.
    int    topModel    = min( scenario, 2);
    double gammat      = (scenario != 3) ? 1.34 : 0.2;
    // Warning:  a small gammatGreen may lead to event weights above unity.
    double gammatGreen = (scenario != 4) ? 1.34 : 0.2;

    // Generator.
    Pythia pythia;

    // Feed in desired values.
    pythia.settings.mode("TopThreshold:model", topModel);
    pythia.particleData.m0(6, mt);
    pythia.readString("6:doForceWidth = true");
    pythia.particleData.mWidth( 6, gammat);
    pythia.settings.parm("TopThreshold:tWidthGreen",   gammatGreen);
    pythia.settings.parm("TopThreshold:thrRegion",     thrRegion);
    pythia.settings.mode("TopThreshold:alphasOrder",   alphasOrder);
    pythia.settings.parm("TopThreshold:alphasValue",   alphasValue);
    pythia.settings.parm("TopThreshold:ggSingletFrac", ggSingletFrac);
    pythia.settings.parm("TopThreshold:qqSingletFrac", qqSingletFrac);

    // Process selection and collision energy.
    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.settings.parm("Beams:eCM", eCM);

    // Restrict to threshold region.
    if (thresholdOnly) {
      pythia.readString("PhaseSpace:mHatMin = 300.");
      pythia.readString("PhaseSpace:mHatMax = 400.");
    }
    else pythia.readString("PhaseSpace:mHatMin = 200.");

    // Switch off (most) code parts not relevant here. Reduce printout.
    if (!fullEvents) {
      pythia.readString("PartonLevel:ISR = off");
      pythia.readString("PartonLevel:FSR = off");
      pythia.readString("PartonLevel:MPI = off");
      pythia.readString("HadronLevel:all = off");
    }
    pythia.readString("Next:numberCount = 100000");

    // If Pythia fails to initialize, abort this scenario.
    // Use critical to avoid thread clashes when initializing in parallel.
    bool ok = true;
#ifdef OPENMP
#pragma omp critical
#endif
    {
      cout << " =============================================================="
        "=========\n Initializing PYTHIA for scenario " << scenario
           << " with topModel " << topModel << " and top mass " << mt
           << " GeV" << endl;
      ok = pythia.init();
      cout << " Completed initialization of PYTHIA for scenario " << scenario
           << "\n ============================================================"
        "===========\n";
    }
    if (!ok) continue;

    // Counters. Real numbers since some events can come with weights > 1.
    double nEventEff = 0., nBelow2mt = 0., nBelow = 0., nBelowgg = 0.,
      nBelowqq = 0., nAbovegg = 0., nAboveqq = 0., eBelow = 0.;

    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events.
      if (!pythia.next()) continue;

      // Original threshold value and other properties.
      double eThr = pythia.info.toponiumE;
      double mHat = pythia.info.mHat();
      double wt   = pythia.info.weight();

      // Histogram threshold value and ttbar invariant mass.
      // Use critical to avoid thread clashes when filling histogram
      // in parallel.
#ifdef OPENMP
#pragma omp critical
#endif
      {
        mHatLow[scenario].fill( mHat, wt);
        mHatLowZoom[scenario].fill( mHat, wt);
        mThrOrig[scenario].fill( eThr, wt);
        mThrOrigZoom[scenario].fill( eThr, wt);
      }

      // Statistics for below- and above-threshold cross sections.
      nEventEff += wt;
      if (mHat < 2. * mt) nBelow2mt += wt;
      if (eThr < 0.) {
        nBelow += wt;
        if (pythia.info.code() == 601) nBelowgg += wt;
        else nBelowqq += wt;
        eBelow += eThr * wt;
      } else {
        if (pythia.info.code() == 601) nAbovegg += wt;
        else nAboveqq += wt;
      }

    // End of event loop.
    }

    // Normalization and statistics for this scenario.
    // Use critical to avoid thread clashes.
#ifdef OPENMP
#pragma omp critical
#endif
    {
      cout << "\n ===== End-of-run statistics for scenario " << scenario
           << " with topModel " << topModel << " and top mass " << mt
           << " GeV =====" << endl;

      // Normalize histogram to cross section, in pb/GeV.
      double sigmaGen         = 1e9 * pythia.info.sigmaGen();
      double sigmapb          = sigmaGen / nEventEff;
      mHatLow[scenario]      *= sigmapb;
      mHatLowZoom[scenario]  *= 4. * sigmapb;
      mThrOrig[scenario]     *= sigmapb;
      mThrOrigZoom[scenario] *= 2.5 * sigmapb;

      // Statistics and cross sections.
      pythia.stat();
      double sigma2mt = sigmapb * nBelow2mt;
      double sigAbv   = sigmapb * (nEventEff - nBelow);
      double sigAbvgg = sigmapb * nAbovegg;
      double sigAbvqq = sigmapb * nAboveqq;
      double sigBel   = sigmapb * nBelow;
      double sigBelgg = sigmapb * nBelowgg;
      double sigBelqq = sigmapb * nBelowqq;
      cout << "\n total sigma in studied region = " << fixed << setprecision(3)
           << sigmaGen << " pb whereof below 2 mt = " << sigma2mt << " pb"
           << endl << " sigma above threshold = "  << sigAbv << " pb" << endl
           << "   whereof gg = " << sigAbvgg << " pb and qq = " << sigAbvqq
           << " pb" << endl << " sigma below threshold = " << sigBel
           << " pb" << endl << "   whereof gg = "
           << sigBelgg << " pb and qq = " << sigBelqq << " pb" << endl;
      eBelow /= max(1., nBelow);
      cout << " average energy for below-threshold part = " << eBelow
           << " GeV" << endl;
    }

  // End loop over threshold cases.
  }

  // No more parallelism beyond this point.
#ifdef OPENMP
#pragma omp barrier
#endif

  // Histograms with pyplot. Common notation.
  HistPlot hpl("plot369");
  string symbol[5] = {"h,black", "h,olive", "h,red", "h,orange", "h,blue"};
  string legend[5] = {"Born", "Coulomb", "Green's, wide BW and ME",
      "Green's, narrow BW", "Green's, narrow ME"};

  // Invariant mass of the ttbar pair.
  hpl.frame("fig369", "invariant mass of ttbar pair, near threshold",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 5; ++scenario)
    hpl.add( mHatLow[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();
  hpl.frame("", "invariant mass of ttbar pair, near threshold, zoom-in",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 5; ++scenario)
    hpl.add( mHatLowZoom[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();

  // The formal E_thr value above/below threshold, key for reweighting.
  hpl.frame("", "Energy above or below threshold",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 5; ++scenario)
    hpl.add( mThrOrig[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();
  hpl.frame("", "Energy above or below threshold, zoom-in",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  for (int scenario = 0; scenario < 5; ++scenario)
    hpl.add( mThrOrigZoom[scenario], symbol[scenario], legend[scenario]);
  hpl.plot();

  // Done.
  return 0;
}
