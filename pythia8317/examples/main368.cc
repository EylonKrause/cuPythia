// main368.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Torbjorn Sjostrand <torbjorn.sjostrand@fysik.lu.se>

// Keywords: top; toponium

// Plot top threshold factors according to chosen parameter settings.
// Based on ttbar Coulomb factors and Green's functions according to FKS:
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613.

#include "Pythia8/Pythia.h"
using namespace Pythia8;

//==========================================================================

// Implementation of the ttbar Coulomb factors and Green's functions in
// V. Fadin,  V. Khoze and T. Sjostrand, Z. Phys. C48 (1990) 613.

class TopThr  {

public:

  // Save input values for given scenario.
  void setup( double mtIn, double gammatIn, int alphasOrder,
    double alphasValue) { mt = mtIn; gammat = gammatIn;
    alphas.init(alphasValue, alphasOrder);}

  // Return threshold factor value, at eNow above or below threshold.
  // isGreen true or false gives Green's function or Coulomb factor.
  // isSinglet gives singlet or octet contribution.
  double value( double eNow, bool isGreen, bool isSinglet) {

    // Negative eNow not allowed in Coulomb.
    if (!isGreen && eNow <= 0.) return 0.;

    // alpha_strong value.
    double q2Thr = pow2(eNow) + pow2(gammat);
    double q2alp = mt * sqrt(q2Thr);
    alps         = alphas.alphaS(q2alp);
    double beta  = sqrtpos(1. - pow2( 2. * mt / (2. * mt + eNow)));

    // Attractive expressions for singlet case.
    if (isSinglet) {
      if (isGreen) return imGreenSin( eNow);
      double xAttr = (4. / 3.) * M_PI * alps / beta;
      return beta * xAttr / (1. - exp(-xAttr));

    // Repulsive expressions for octet case.
    } else {
      if (isGreen) return imGreenOct( eNow);
      double xRepu = (1. / 6.) * M_PI * alps / beta;
      return beta * xRepu / (exp(xRepu) - 1);
    }

  }

  // Imaginary part of Green's function for singlet state.
  // Factor m_t^2 / (4 pi) omitted since cancelled when applied for sigma.
  double imGreenSin(double eNow) {

    // Basic expressions.
    double ps   = (2. / 3.) * mt * alps;
    double egrt = sqrt(eNow * eNow + gammat * gammat);
    double p1   = sqrt( 0.5 * mt * (egrt - eNow));
    double p2   = sqrt( 0.5 * mt * (egrt + eNow));

    // Sum over resonance contributions.
    double ressum = 0.;
    for (int n = 1; n < 21; ++n)
      ressum += (gammat * ps * n + p2 * (n*n * egrt + ps * ps / mt))
      / ( pow4(n) * (pow2(eNow + ps * ps / (mt * n*n)) + gammat * gammat) );

    // Combine with non-resonant terms and done.
    return p2 / mt + (2. * ps / mt) * atan(p2 / p1)
      + 2. * pow2(ps / mt) * ressum;
  }

  // Imaginary part of Green's function for octet state.
  // Factor m_t^2 / (4 pi) omitted since cancelled when applied for sigma.
  double imGreenOct(double eNow) {

    // Basic expressions.
    double p8   = - (1. / 12.) * mt * alps;
    double egrt = sqrt(eNow * eNow + gammat * gammat);
    double p1   = sqrt( 0.5 * mt * (egrt - eNow));
    double p2   = sqrt( 0.5 * mt * (egrt + eNow));

    // Sum over resonance contributions.
    double ressum = 0.;
    for (int n = 1; n < 21; ++n)
      ressum += mt * p2 / (pow2(n * p1 - p8) + pow2(n * p2));

    // Combine with non-resonant terms and done.
    return p2 / mt + (2. * p8 / mt) * atan(p2 / p1)
      + 2. * pow2(p8 / mt) * ressum;
  }

 private:

  // Commonly available variables.
  double mt, gammat, alps;

  // Need alphaStrong with special scale.
  AlphaStrong alphas;

};

//==========================================================================

int main() {

  // useOld = true to reproduce FKS figures, else use more modern values.
  bool useOld = false;

  // showMod = true shows Green's function also with modified argument.
  // May give too busy figures, so optional.
  bool showMod = false;

  // splitFigs = true gives each figure separately, else as one single file.
  bool splitFigs = false;

  // Generator.
  Pythia pythia;

  // Parameters to set.
  int    alphasOrder;
  double mt, gammat, alphasValue;

  // Old setup from FKS or a more relevent current setup.
  if (useOld) {
    mt            = 200.;
    gammat        = 2.5;
    alphasOrder   = 1;
    alphasValue   = 0.134;
  } else {
    mt            = 172.5;
    gammat        = 1.34;
    alphasOrder   = 2;
    alphasValue   = 0.118;
  }

  // Create and initialize top threshold object.
  TopThr topThr;
  topThr.setup( mt, gammat, alphasOrder, alphasValue);

  // Histograms, for narrower or broader E range around threshold.
  Hist greensin("imGreen singlet", 200, -10., 30.);
  Hist greenoct("imGreen octet", 200, -10., 30.);
  Hist greensinMod("imGreen singlet", 200, -10., 30.);
  Hist greenoctMod("imGreen octet", 200, -10., 30.);
  Hist betaval("beta threshold factor",200, -10., 30.);
  Hist coulsin("singlet Coulomb factor",200, -10., 30.);
  Hist couloct("octet Coulomb factor",200, -10., 30.);
  Hist greensin2("imGreen singlet", 240, -20., 100.);
  Hist greenoct2("imGreen octet", 240, -20., 100.);
  Hist greensin2Mod("imGreen singlet", 240, -20., 100.);
  Hist greenoct2Mod("imGreen octet", 240, -20., 100.);
  Hist betaval2("beta threshold factor",240, -20., 100.);
  Hist coulsin2("singlet Coulomb factor",240, -20., 100.);
  Hist couloct2("octet Coulomb factor",240, -20., 100.);

  // Loop over energies to plot - narrower range.
  for (int iE = 0; iE < 200; ++iE) {
    double eNow  = -9.9 + 0.2 * iE;
    double beta2 = 1. - pow2( 2. * mt / (2. * mt + eNow));

    // Green's function expressions.
    double valNowS = topThr.value( eNow, true, true);
    greensin.fill( eNow, valNowS);
    double valNowO = topThr.value( eNow, true, false);
    greenoct.fill( eNow, valNowO);

    // Relativistic argument to Green's function is equivalent at threshold.
    double eMod  = mt * beta2;
    double valNowSMod = topThr.value( eMod, true, true);
    greensinMod.fill( eNow, valNowSMod);
    double valNowOMod = topThr.value( eMod, true, false);
    greenoctMod.fill( eNow, valNowOMod);

    // For positive energies also beta and Coulomb.
    if (eNow > 0.) {
      double beta = sqrt(beta2);
      betaval.fill( eNow, beta);
      double fAttr = topThr.value( eNow, false, true);
      coulsin.fill( eNow, fAttr);
      double fRepu = topThr.value( eNow, false, false);
      couloct.fill( eNow, fRepu);
    }
  }

  // Loop over energies to plot - broader range.
  for (int iE = 0; iE < 240; ++iE) {
    double eNow = -19.75 + 0.5 * iE;
    double beta2 = 1. - pow2( 2. * mt / (2. * mt + eNow));

    // Green's function expressions.
    double valNowS = topThr.value( eNow, true, true);
    greensin2.fill( eNow, valNowS);
    double valNowO = topThr.value( eNow, true, false);
    greenoct2.fill( eNow, valNowO);

    // Relativistic argument to Green's function is equivalent at threshold.
    double eMod  = mt * beta2;
    double valNowSMod = topThr.value( eMod, true, true);
    greensin2Mod.fill( eNow, valNowSMod);
    double valNowOMod = topThr.value( eMod, true, false);
    greenoct2Mod.fill( eNow, valNowOMod);

    // For positive energies also beta and Coulomb.
    if (eNow > 0.) {
      double beta = sqrt(beta2);
      betaval2.fill( eNow, beta);
      double fAttr = topThr.value( eNow, false, true);
      coulsin2.fill( eNow, fAttr);
      double fRepu = topThr.value( eNow, false, false);
      couloct2.fill( eNow, fRepu);
    }
  }

  // Histograms with pyplot. Frame dimensions.
  HistPlot hpl("plot368");
  double width  = 4.8;
  double height = 4.8;

  // Singet and octet contribution in -10 < E < 30.
  hpl.frame(splitFigs ? "fig368narrow" : "fig368", "", "$E$ (GeV)",
    "rate (arbitrary units)", width, height);
  hpl.add( betaval, "-,black", "pure beta threshold");
  hpl.add( coulsin, "--,blue", "singlet Coulomb factor");
  hpl.add( couloct, "--,cyan", "octet Coulomb factor");
  hpl.add( greensin, "--,red", "singlet Green's function");
  hpl.add( greenoct, "--,magenta", "octet Green's function");
  if (showMod) {
    hpl.add( greensinMod, "-.,red", "modified singlet Green's function");
    hpl.add( greenoctMod, "-.,magenta", "modified octet Green's function");
  }
  hpl.plot();

  // Singet and octet contribution in -20 < E < 100.
  hpl.frame(splitFigs ? "fig368wide" : "", "", "$E$ (GeV)",
    "rate (arbitrary units)", width, height);
  hpl.add( betaval2, "-,black", "pure beta threshold");
  hpl.add( coulsin2, "--,blue", "singlet Coulomb factor");
  hpl.add( couloct2, "--,cyan", "octet Coulomb factor");
  hpl.add( greensin2, "--,red", "singlet Green's function");
  hpl.add( greenoct2, "--,magenta", "octet Green's function");
  if (showMod) {
    hpl.add( greensin2Mod, "-.,red", "modified singlet Green's function");
    hpl.add( greenoct2Mod, "-.,magenta", "modified octet Green's function");
  }
  hpl.plot();

  // Singlet contribution only, for direct comparison with FKS figure.
  hpl.frame(splitFigs ? "fig368singlet" : "", "",  "$E$ (GeV)",
    "rate (arbitrary units)", width, height);
  hpl.add( betaval, "-,black", "pure beta threshold");
  hpl.add( coulsin, "--,blue", "singlet Coulomb factor");
  hpl.add( greensin, "-.,red", "singlet Green's function");
  if (showMod && !useOld) hpl.add( greensinMod, "-.,blue",
    "modified singlet Green's function");
  if (useOld) hpl.plot(-10., 30., 0., 0.74);
  else hpl.plot();

  // Octet contribution only, for direct comparison with FKS figure.
  hpl.frame(splitFigs ? "fig368octet" : "", "",  "$E$ (GeV)",
    "rate (arbitrary units)", width, height);
  hpl.add( betaval, "-,black", "pure beta threshold");
  hpl.add( couloct, "--,cyan", "octet Coulomb factor");
  hpl.add( greenoct, "-.,magenta", "octet Green's function");
  if (showMod && !useOld) hpl.add( greenoctMod, "-.,blue",
    "octet Green's function");
  if (useOld) hpl.plot(-10., 30., 0., 0.37);
  else hpl.plot();

  // Done.
  return 0;
}

//============================================================================
