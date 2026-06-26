// KinematicVariations.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Stephen Mrenna, Christian Bierlich, Philip Ilten,
// Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#ifndef Pythia8_KinematicVariations_H
#define Pythia8_KinematicVariations_H

// Include Pythia headers.
#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// Calculate the weight for an event given a different set of flavor
// parameters used in the hadronization.

class KinematicVariations {

public:

  // Constructor, given an intialized Pythia object.
  KinematicVariations(Settings &settings) : KinematicVariations(
    settings.parm("StringZ:aLund"),
    settings.parm("StringZ:bLund"),
    settings.parm("StringZ:rFactC"),
    settings.parm("StringZ:rFactB"),
    settings.parm("StringPT:sigma"),
    settings.parm("VariationFrag:zHead")) {}

  // Constructor, given the default base parameters.
  KinematicVariations(double aLund, double bLund, double rFactC,
    double rFactB, double sigma, double zHead) :
    pythia("", false) {
    pythia.settings.flag("ProcessLevel:all", false);
    pythia.settings.flag("Print:quiet", true);
    pythia.settings.parm("StringZ:aLund", aLund);
    pythia.settings.parm("StringZ:bLund", bLund);
    pythia.settings.parm("StringZ:rFactC", rFactC);
    pythia.settings.parm("StringZ:rFactB", rFactB);
    pythia.settings.parm("StringPT:sigma", sigma);
    pythia.settings.parm("VariationFrag:zHead", zHead);
    pythia.init();
  }

  // Write string breaks.
  string write(const vector<int>& breaks) {
    string out = "{";
    for (const int& val : breaks) out += toString(val) + ",";
    return out.substr(0, out.length() - 1) + "}";}
  string write(const vector<double>& breaks) {
    string out = "{";
    for (const double& val : breaks) out += toString(val) + ",";
    return out.substr(0, out.length() - 1) + "}";}

  // Calculate the full kinematic weight.
  double weight(double aLund, double bLund, double rFactC, double rFactB,
    double sigma, const vector<int>& zIntBreaks,
    const vector<double>& zDblBreaks, const vector<double>& pTBreaks) {
    return weight(aLund, bLund, rFactC, rFactB, zIntBreaks, zDblBreaks) *
      weight(sigma, pTBreaks);
  }

  // Calculate just the z weight.
  double weight(double aLund, double bLund, double rFactC, double rFactB,
    const vector<int>& zIntBreaks, const vector<double>& zDblBreaks) {
    double wgt = 1.;
    for (int iBrk = 0; iBrk < int(zIntBreaks.size()/2); ++iBrk) {
      wgt *= pythia.info.weightContainerPtr
        ->weightsFragmentation.zWeight(aLund, bLund, rFactC, rFactB,
          zIntBreaks[2*iBrk + 0], zIntBreaks[2*iBrk + 1],
          zDblBreaks[3*iBrk + 0], zDblBreaks[3*iBrk + 1],
          zDblBreaks[3*iBrk + 2]);
    }
    return wgt;
  }

  // Calculate just the pT weight.
  double weight(double sigma, const vector<double>& pTBreaks) {
    double wgt = 1.;
    for (int iBrk = 0; iBrk < int(pTBreaks.size()/2); ++iBrk)
      wgt *= pythia.info.weightContainerPtr
        ->weightsFragmentation.pTWeight(
          sigma, pTBreaks[2*iBrk + 0], pTBreaks[2*iBrk + 1]);
    return wgt;
  }

private:

  // Pythia object.
  Pythia pythia;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_KinematicVariations_H
