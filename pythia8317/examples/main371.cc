// main371.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: top; toponium

// Study ttbar production above and below threshold,
// starting out from the equations in
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613,
// with some further assumptions in the below-threshold simulation.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//============================================================================

int main() {

  // Restrict to threshold region only or not.
  bool thresholdOnly = true;

  // Full event generation or only cross-section-oriented studies.
  bool fullEvents = false;

  // Number of events.
  int nEvent = 1000000;

  // LHC collision energy.
  double eCM = 13000.;

  // Example of possible setup values.
  // topModelSetup = 0 is Born, = 1 is simple Coulomb, = 2 is Green's,
  // = 3 is reduced top quark width, = 4 reduced width in cross section.
  int    topModelSetup   = 3;

  // Suitable values for given topModelSetup.
  int    topModel        = min( topModelSetup, 2);
  double mt              = 172.5;
  double gammat          = (topModelSetup != 3) ? 1.34 : 0.1;
  // Warning:  a small gammatGreen may lead to event weights above unity.
  double gammatGreen     = (topModelSetup != 4) ? 1.34 : 0.2;
  double thresholdRegion = 10.;
  int    alphasOrder     = 2;
  double alphasValue     = 0.118;
  double ggSingletFrac   = 2./7.;
  double qqSingletFrac   = 0.;
  bool   pseudoscalar    = true;
  double eBegDamp        = 10.;
  double eEndDamp        = 20.;

  // Generator.
  Pythia pythia;

  // Process selection. LHC at 13 TeV.
  pythia.readString("Top:gg2ttbar = on");
  pythia.readString("Top:qqbar2ttbar = on");
  pythia.settings.parm("Beams:eCM", eCM);

  // Restrict to threshold region or not.
  if (thresholdOnly) {
    pythia.readString("PhaseSpace:mHatMin = 300.");
    pythia.readString("PhaseSpace:mHatMax = 400.");
  }
  else pythia.readString("PhaseSpace:mHatMin = 200.");

  // Feed in desired values.
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
  pythia.settings.flag("TopThreshold:pseudoscalar", pseudoscalar);
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

  // Histograms.
  Hist mHatOri("original threshold", 100, -20., 80.);
  Hist mHatOriLow("original threshold", 100, -10., 30.);
  Hist mHatThr("corrected threshold", 100, -20., 80.);
  Hist mHatLow("ttbar invariant mass, low", 100, 300., 400.);
  Hist mHatAll("ttbar invariant mass, all", 100, 200., 1200.);
  Hist mTopMax("max(m_t, m_tbar) mass distribution", 100, 100., 200.);
  Hist mTopMin("min(m_t, m_tbar) mass distribution", 100, 100., 200.);
  Hist betaPair("beta factor of pair", 100, 0., 1.);
  Hist mHatOriBT("original threshold, below", 100, -20., 80.);
  Hist mHatOriBTLow("original threshold, below", 100, -10., 30.);
  Hist mHatThrBT("corrected threshold, below", 100, -10., 30.);
  Hist mHatLowBT("ttbar invariant mass, low, below", 100, 300., 400.);
  Hist mHatAllBT("ttbar invariant mass, all, below", 100, 200., 1200.);
  Hist mTopMaxBT("max(m_t, m_tbar) mass distribution, below", 100, 100., 200.);
  Hist mTopMinBT("min(m_t, m_tbar) mass distribution, below", 100, 100., 200.);
  Hist betaPairBT("beta factor of pair, below", 100, 0., 1.);
  Hist mShiftBT("shift in top mass, after - before", 100, -45., 5.);

  // Counters. Real numbers since some events can come with weights > 1.
  double nEventEff = 0., nBelow2mt = 0., nBelow = 0., nBelowgg = 0.,
    nBelowqq = 0., nAbovegg = 0., nAboveqq = 0., eBelow = 0.;

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events.
    if (!pythia.next()) continue;
    double wt   = pythia.info.weight();

    // Original threshold value and other properties.
    double eThr = pythia.info.toponiumE;
    double mt1T = pythia.info.toponiumm3;
    double mt2T = pythia.info.toponiumm4;
    double mHat = pythia.info.mHat();
    double mt1  = pythia.info.m3Hat();
    double mt2  = pythia.info.m4Hat();
    double rat1 = pow2( mt1 / mHat);
    double rat2 = pow2( mt2 / mHat);
    double beta = sqrtpos( pow2(1 - rat1 - rat2) - 4. * rat1 * rat2);

    // Histogram for all events.
    mHatOri.fill( eThr, wt);
    mHatOriLow.fill( eThr, wt);
    mHatThr.fill( mHat - mt1 - mt2, wt);
    mHatLow.fill( mHat, wt);
    mHatAll.fill( mHat, wt);
    mTopMax.fill( max(mt1,mt2), wt);
    mTopMin.fill( min(mt1,mt2), wt);
    betaPair.fill( beta, wt);

    // Histogram for below-threshold events (and some above).
    nEventEff += wt;
    if (mHat < 2. * mt) nBelow2mt += wt;
    if (eThr < 0.) {
      nBelow += wt;
      if (pythia.info.code() == 601) nBelowgg += wt;
      else nBelowqq += wt;
      eBelow += eThr * wt;
      mHatOriBT.fill( eThr, wt);
      mHatOriBTLow.fill( eThr, wt);
      mHatThrBT.fill( mHat - mt1 - mt2, wt);
      mHatLowBT.fill( mHat, wt);
      mHatAllBT.fill( mHat, wt);
      mTopMaxBT.fill( max(mt1,mt2), wt);
      mTopMinBT.fill( min(mt1,mt2), wt);
      betaPairBT.fill( beta, wt);
      mShiftBT.fill( mt1 - mt1T, wt);
      mShiftBT.fill( mt2 - mt2T, wt);
    } else {
      if (pythia.info.code() == 601) nAbovegg += wt;
      else nAboveqq += wt;
    }

  // End of event loop.
  }

  // Normalize histogram to cross section, in pb/GeV.
  double sigmaGen = 1e9 * pythia.info.sigmaGen();
  double sigmapb = sigmaGen  / nEventEff;
  mHatOri      *= sigmapb;
  mHatOriLow   *= 2.5 * sigmapb;
  mHatThr      *= sigmapb;
  mHatLow      *= sigmapb;
  mHatAll      *= 0.1 * sigmapb;
  mTopMax      *= sigmapb;
  mTopMin      *= sigmapb;
  betaPair     *= 100. * sigmapb;
  mHatOriBT    *= sigmapb;
  mHatOriBTLow *= 2.5 * sigmapb;
  mHatThrBT    *= sigmapb;
  mHatLowBT    *= sigmapb;
  mHatAllBT    *= 0.1 * sigmapb;
  mTopMaxBT    *= sigmapb;
  mTopMinBT    *= sigmapb;
  betaPairBT   *= 100. * sigmapb;
  mShiftBT     *= 2. * sigmapb;

  // Statistics and cross sections
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
       << " pb" << endl << "   whereof gg = " << sigBelgg << " pb and qq = "
       << sigBelqq << " pb" << endl;
  eBelow /= max(1., nBelow);
  cout << " average energy for below-threshold part = " << eBelow
       << " GeV" << endl;

  // Histograms with pyplot.
  HistPlot hpl("plot371");
  bool doOptional = true;

  // Spectra relative to event-by-event threshold, including below.
  hpl.frame("fig371", "Energy above/below threshold, original tops",
    "$E$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}E$ (pb/GeV)");
  hpl.add( mHatOri - mHatOriBT, "h,black", "above-threshold part");
  hpl.add( mHatOriBT, "h,red", "below-threshold part");
  hpl.plot();

  // Ditto, but with new t/tbar masses for below-threshold part.
  hpl.frame("", "Energy above threshold, new top masses where necessary",
    "$E'$ (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}E'$ (pb/GeV)");
  hpl.add( mHatThr, "h,black", "all");
  hpl.add( mHatThrBT, "h,red", "below-threshold part");
  hpl.plot();

  // Pair mass spectra close to threshold.
  hpl.frame("", "Invariant mass of ttbar pair, near threshold",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  hpl.add( mHatLow, "h,black", "all");
  hpl.add( mHatLowBT, "h,red", "below-threshold part");
  hpl.plot();

  // Pair mass spectra over full range.
  hpl.frame("", "Invariant mass of ttbar pair, full range",
    "$m(t+tbar) (GeV)$", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)");
  hpl.add( mHatAll, "h,black", "all");
  hpl.add( mHatAllBT, "h,red", "below-threshold part");
  hpl.plot();

  // Top/antitop mass spectra.
  hpl.frame("", "Larger of top and antitop masses",
    "max($m(t), m(tbar)$) (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ "
    "(pb/GeV)");
  hpl.add( mTopMax, "h,black", "all");
  hpl.add( mTopMaxBT, "h,red", "below-threshold part");
  hpl.plot(100., 200., 0.01, 100., true, false);

  hpl.frame("", "Smaller of top and antitop masses",
    "min($m(t), m(tbar)$) (GeV)", "$\\mathrm{d}\\sigma/\\mathrm{d}m$ "
    "(pb/GeV)");
  hpl.add( mTopMin, "h,black", "all");
  hpl.add( mTopMinBT, "h,red", "below-threshold part");
  hpl.plot(100., 200., 0.01, 100., true, false);

  // Top/antitop separation.
  hpl.frame("", "Beta factor in t-tbar pair",
    "$\\beta$", "$\\mathrm{d}\\sigma/\\mathrm{d}\\beta$ (pb)");
  hpl.add( betaPair, "h,black", "all");
  hpl.add( betaPairBT, "h,red", "below-threshold part");
  hpl.plot();

  // Redefinition change of top masses in below-threshold pairs.
  hpl.frame("", "Shift of below-threshold top masses", "$\\Delta m$ (GeV)",
    "$\\mathrm{d}\\sigma/\\mathrm{d}(\\Delta m)$ (pb/GeV)");
  hpl.add( mShiftBT, "h,red", "below-threshold part");
  hpl.plot();

  // For inclusion in talk.
  if (doOptional) {
    hpl.frame("fig371mass", "", "$\\hat{m}$ (GeV)",
      "$\\mathrm{d}\\sigma/\\mathrm{d}\\hat{m}$ (pb/GeV)", 4.8, 4.8);
    hpl.add( mHatLowBT, "h,red", "below threshold");
    hpl.add( mHatLow, "h,black", "all");
    hpl.plot();
    hpl.frame("fig371BreitWigners", "", "$m$ (GeV)",
      "$\\mathrm{d}\\sigma/\\mathrm{d}m$ (pb/GeV)", 4.8, 4.8);
    hpl.add( mTopMax, "h,black", "max($m_{t1},m_{t2}$) all");
    hpl.add( mTopMin, "h,blue", "min($m_{t1},m_{t2}$) all");
    hpl.add( mTopMaxBT, "h,red", "max($m_{t1},m_{t2}$) below threshold");
    hpl.add( mTopMinBT, "h,magenta", "min($m_{t1},m_{t2}$) below threshold");
    hpl.plot(140., 200., 0.01, 1000., true, false);
  }

  // Done.
  return 0;
}
