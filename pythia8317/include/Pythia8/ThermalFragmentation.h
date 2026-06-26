// ThermalFragmerntation.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains helper classes for thermal fragmentation.
// ThermalStringFlav is used to select quark and hadron flavours.
// ThermalStringPT is used to select transverse momenta.
// ThermalStringEnd described fragmentation from one end of the string.
// ThermalFragmentation is the top-level class for this model.

#ifndef Pythia8_ThermalFragmentation_H
#define Pythia8_ThermalFragmentation_H

#include "Pythia8/FragmentationModel.h"
#include "Pythia8/MiniStringFragmentation.h"
#include "Pythia8/StringFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The ThermalStringFlav class is used to select quark and hadron flavours.

class ThermalStringFlav : public StringFlav {

public:

  // Constructor.
  ThermalStringFlav() : mesonNonetL1(), temperature(), tempPreFactor(),
    nNewQuark(), mesMixRate1(), mesMixRate2(), mesMixRate3(),
    baryonOctWeight(), baryonDecWeight(), hadronConstIDs(), possibleHadrons(),
    possibleRatePrefacs(), possibleHadronsLast(), possibleRatePrefacsLast(),
    hadronIDwin(), quarkIDwin(), hadronMassWin() {}

  // Destructor.
  ~ThermalStringFlav() {}

  // Initialize data members.
  void init() override;

  // Pick a new flavour (including diquarks) given an incoming one.
  FlavContainer pick(FlavContainer& flavOld,
    double pT, double kappaModifier, bool allowPop) override;

  // Return chosen hadron in case of thermal model.
  virtual int getHadronIDwin() { return hadronIDwin; }

  // Return hadron mass. Used one if present, pick otherwise.
  double getHadronMassWin(int idHad) override { return
    ((hadronMassWin < 0.0) ? particleDataPtr->mSel(idHad) : hadronMassWin); }

  // Combine two flavours into hadron for last two remaining flavours
  // for thermal model.
  int combineLastThermal(FlavContainer& flav1, FlavContainer& flav2,
    double pT, double kappaModifier);

  // Return already set hadron id or combination of the two flavours.
  int getHadronID(FlavContainer& flav1, FlavContainer& flav2, double pT = -1.0,
    double kappaModifier = -1.0, bool finalTwo = false) override {
    if (finalTwo) return combineLastThermal(flav1, flav2, pT, kappaModifier);
    if ( (hadronIDwin != 0) && (quarkIDwin != 0)) return getHadronIDwin();
    return combine(flav1, flav2); }

  // Check if quark-diquark combination should be added. If so add.
  void addQuarkDiquark(vector< pair<int,int> >& quarkCombis,
    int qID, int diqID, int hadronID) {
    bool allowed = true;
    for (int iCombi = 0; iCombi < int(quarkCombis.size()); iCombi++)
      if ( (qID   == quarkCombis[iCombi].first ) &&
           (diqID == quarkCombis[iCombi].second) ) allowed = false;
    if (allowed) quarkCombis.push_back( (hadronID > 0) ?
      make_pair( qID,  diqID) : make_pair(-qID, -diqID) ); }

protected:

  // Settings for thermal model.
  bool   mesonNonetL1;
  double temperature, tempPreFactor;
  int    nNewQuark;
  double mesMixRate1[2][6], mesMixRate2[2][6], mesMixRate3[2][6];
  double baryonOctWeight[6][6][6][2], baryonDecWeight[6][6][6][2];

  // Key = hadron id, value = list of constituent ids.
  map< int, vector< pair<int,int> > > hadronConstIDs;
  // Key = initial (di)quark id, value = list of possible hadron ids
  //                                     + nr in hadronConstIDs.
  map< int, vector< pair<int,int> > > possibleHadrons;
  // Key = initial (di)quark id, value = prefactor to multiply rate.
  map< int, vector<double> > possibleRatePrefacs;
  // Similar, but for combining the last two (di)quarks. Key = (di)quark pair.
  map< pair<int,int>, vector< pair<int,int> > > possibleHadronsLast;
  map< pair<int,int>, vector<double> > possibleRatePrefacsLast;

  // Selection in thermal model.
  int    hadronIDwin, quarkIDwin;
  double hadronMassWin;

};

//==========================================================================

// The ThermalStringPT class is used to select transverse momenta.

class ThermalStringPT : public StringPT {

public:

  // Constructor.
  ThermalStringPT() : temperature(), tempPreFactor(), fracSmallX() {}

  // Destructor.
  ~ThermalStringPT() {}

  // Initialize data members.
  void init() override;

  // Return px and py as a pair.
  pair<double, double> pxy(int idIn = 0, double kappaModifier = -1.0) override;

  // Exponential suppression of given pT2; used in MiniStringFragmentation.
  double suppressPT2(double pT2) override {
    return exp(-sqrt(pT2)/temperature); }

protected:

  // Initialization data, to be read from Settings.
  double temperature, tempPreFactor, fracSmallX;

private:

  // Evaluate Bessel function K_{1/4}(x).
  double BesselK14(double x);

};

//==========================================================================

// The ThermalFragmentation class handles the alternative fragmentation,
// using both the StringFragmentation and MiniStringFragmentation classes.

class ThermalFragmentation : public FragmentationModel {

public:

  // Constructor (creates string, string-end and mini-string pointers).
  ThermalFragmentation();

  // Destructor (deletes string, string-end and mini-string pointers).
  ~ThermalFragmentation() override;

  // Initialize and save pointers.
  bool init(StringFlav* flavSelPtrIn = nullptr,
    StringPT* pTSelPtrIn = nullptr, StringZ* zSelPtrIn = nullptr,
    FragModPtr fragModPtrIn = nullptr) override;

  // Do the fragmentation: driver routine.
  bool fragment(int iSub, ColConfig& colConfig, Event& event,
    bool isDiff = false, bool systemRecoil = true) override;

  // Classes for flavour, pT and z generation.
  ThermalStringFlav* thermalFlavSelPtr{};
  ThermalStringPT*   thermalPTSelPtr{};
  StringZ*           zSelPtr{};

  // Internal StringFragmentation and MiniStringFragmentation objects.
  StringFragmentation* stringFragPtr{};
  MiniStringFragmentation* ministringFragPtr{};

private:

  // Parameters controlling the fragmentation.
  double mStringMin{};
  bool tryMiniAfterFailedFrag{};

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_ThermalFragmentation_H
