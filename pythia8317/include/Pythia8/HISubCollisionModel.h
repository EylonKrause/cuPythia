// HISubCollisionModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the ImpactParameterGenerator,
// SubCollision, and SubCollisionModel classes, as well as a set of
// subclasses of SubCollisionModel.
//
// ImpactParameterGenerator: distributes nuclei in impact parameter space.
// SubCollision: a collision between a projectile and a target Nucleon.
// SubCollisionModel: Models the collision probabilities of nucleons.
// BlackSubCollisionModel: A very simple SubCollisionModel.
// NaiveSubCollisionModel: A very simple SubCollisionModel.
// DoubleStrikmanSubCollisionModel: A more advanced SubCollisionModel.

#ifndef Pythia8_HISubCollisionModel_H
#define Pythia8_HISubCollisionModel_H

#include "Pythia8/Pythia.h"
#include "Pythia8/HINucleusModel.h"

namespace Pythia8 {

//==========================================================================

// SubCollision represents a possible collision between a projectile
// and a target nucleon.

class SubCollision {

public:

  // This defines the type of a binary nucleon collision.
  enum CollisionType {
    NONE,       // This is not a collision.
    ELASTIC,    // This is an elastic scattering
    SDEP,       // The projectile is diffractively excited.
    SDET,       // The target is diffractively excited.
    DDE,        // Both projectile and target are diffractively excited.
    CDE,        // Both excited but with central diffraction.
    ABS,        // This is an absorptive (non-diffractive) collision.
    LEXC,       // (Low energy) excitation of target and projectile.
    LANN,       // (Low energy) annihilation.
    LRES        // (Low energy) target--projectile resonance.
  };

  // Constructor with configuration.
  SubCollision(Nucleon & projIn, Nucleon & targIn,
               double bIn, double bpIn, double olappIn, CollisionType typeIn)
    : proj(&projIn), targ(&targIn), b(bIn), bp(bpIn), olapp(olappIn),
    type(typeIn), failed(false) {}

  // Default constructor.
  SubCollision()
    : proj(0), targ(0), b(0.0), bp(0.0), olapp(0.0),
      type(NONE), failed(false) {}

  // Used to order sub-collisions in a set.
  bool operator< (const SubCollision& s) const { return b < s.b; }

  // Return 0 if neither proj or target are neutrons, 1 if target is
  // neutron, 2 if projectile is neutron, and 3 if both are neutrons.
  int nucleons() const {return ( abs(targ->id()) == 2112? 1: 0 ) +
      ( abs(proj->id()) == 2112? 2: 0 );}

  // The projectile nucleon.
  Nucleon* proj;

  // The target nucleon.
  Nucleon* targ;

  // The impact parameter distance between the nucleons in femtometer.
  double b;

  // The impact parameter distance between the nucleons scaled like
  // in Pythia to have unit average for non-diffractive collisions.
  double bp;

  // The overlap for the given impact parameter scaled like in Pythia
  // to have unit average for non-diffractive collisions.
  double olapp;

  // The type of collision.
  CollisionType type;

  // Whether the subcollision failed, i.e. has a failed excitation.
  mutable bool failed;

};

//==========================================================================

// The SubCollisionSet gives a set of subcollisions between two nuclei.

class SubCollisionSet {

public:

  // Default constructor.
  SubCollisionSet() = default;

  // Constructor with subcollisions.
  SubCollisionSet(multiset<SubCollision> subCollisionsIn, double TIn)
    : subCollisionsSave(subCollisionsIn), TSave({TIn, TIn, TIn, TIn}) {}

  // Constructor with subcollisions. Special case for subcollision
  // models with 2x2 fluctuating radii.
  SubCollisionSet(multiset<SubCollision> subCollisionsIn, double TIn,
                  double T12In, double T21In, double T22In)
    : subCollisionsSave(subCollisionsIn), TSave({TIn, T12In, T21In, T22In}) {}

  // Reset the subcollisions.
  bool empty() const { return subCollisionsSave.empty(); }

  // The full elastic amplitude, optionally returning alternate states
  // to gauge fluctuations.
  double T(unsigned i = 0) const { return TSave[i]; }

  // Iterators over the subcollisions.
  multiset<SubCollision>::const_iterator begin() const {
    return subCollisionsSave.begin(); }
  multiset<SubCollision>::const_iterator end() const {
    return subCollisionsSave.end(); }

private:

  // Saved subcollisions.
  multiset<SubCollision> subCollisionsSave;

  // The full elastic amplitude together with alternate states gauging
  // fluctuations.
  vector<double> TSave = {};

};

//==========================================================================

// The SubCollisionModel is is able to model the collision between two
// nucleons to tell which type of collision has occurred. The model
// may manipulate the corresponding state of the nucleons.

class SubCollisionModel {

public:

  // Internal class to report cross section estimates.
  struct SigEst {
    // The cross sections (tot, nd, dd, sdp, sdt, cd, el, bslope).
    vector<double> sig;

    // The estimated error (squared)
    vector<double> dsig2;

    // Which cross sections were actually fitted
    vector<bool> fsig;

    // The estimate of the average (and squared error) impact
    // parameter for inelastic non-diffractive collisions.
    double avNDb, davNDb2, avNOF, davNOF2;

    // Constructor for zeros.
    SigEst(): sig(8, 0.0), dsig2(8, 0.0), fsig(8, false),
              avNDb(0.0), davNDb2(0.0), avNOF(0.0), davNOF2(0.0) {}

  };

  // The default constructor is empty.
  // The avNDb has units femtometer.
  SubCollisionModel(int nParm)
    : sigTarg(11, 0.0), sigErr(8, 0.05), sigTargNN(4, vector<double>(11, 0.0)),
      parmSave(nParm), NInt(100000), NPop(20), sigFuzz(0.2), impactFudge(1),
      fitPrint(true), eCMlow(20.0), avNDb(1.0), avNDolap(1.0),
      projPtr(), targPtr(), sigTotPtr(), sigCmbPtr(), settingsPtr(),
      infoPtr(), rndmPtr(), loggerPtr(), particleDataPtr(),
      idASave(0), idBSave(0), doVarECM(false), doVarBeams(false),
      eMin(0.0), eMax(0.0), eSave(0.0), eCMPts(5), idAList(),
      subCollParmsPtr(), subCollParmsMap(), elasticMode(1),
      elasticFudge(0.0) {}

  // Virtual destructor.
  virtual ~SubCollisionModel() {}

  // Create a new SubCollisionModel of the given model.
  static shared_ptr<SubCollisionModel> create(int model);

  // Virtual init method.
  virtual bool init(int idAIn, int idBIn, double eCMIn);

  // Initialize the pointers.
  void initPtr(NucleusModel & projIn, NucleusModel & targIn,
               SigmaTotal & sigTotIn, Settings & settingsIn,
               Info & infoIn, Rndm & rndmIn) {
    projPtr = &projIn;
    targPtr = &targIn;
    sigTotPtr = &sigTotIn;
    settingsPtr = &settingsIn;
    infoPtr = &infoIn;
    rndmPtr = &rndmIn;
    loggerPtr = infoIn.loggerPtr;
    particleDataPtr = infoIn.particleDataPtr;
  }

  // Initialize low energy treatment.
  void initLowEnergy(SigmaCombined * sigmaCombPtrIn) {
    lowEnergyCache.sigCmbPtr = sigCmbPtr = sigmaCombPtrIn;
    lowEnergyCache.particleDataPtr = particleDataPtr;
  }

  bool hasXSec() const {
    if ( elasticMode )
      return ( sigTargNN[0][0] + sigTargNN[1][0] +
               sigTargNN[2][0] + sigTargNN[3][0]  > 0.0 );
    for ( int i = 0; i< 4; ++i )
      for ( int j = 1; j < 11; ++j ) {
        if ( j == 6 || j == 7 ) continue;
        if ( sigTargNN[i][j] > 0.0 ) return true;
      }
    return false;
  }

  // Access the nucleon-nucleon cross sections assumed
  // for this model.

  // The target total nucleon-nucleon cross section.
  double sigTot() const { return sigTarg[0]; }

  // The target elastic cross section.
  double sigEl() const { return sigTarg[6]; }

  // The target central diffractive excitation cross section.
  double sigCDE() const { return sigTarg[5]; }

  // The target single diffractive excitation cross section (both sides).
  double sigSDE() const { return sigSDEP() + sigSDET(); }

  // The target single diffractive excitation cross section (projectile).
  double sigSDEP() const { return sigTarg[3] - sigTarg[1] - sigTarg[2]; }

  // The target single diffractive excitation cross section (target).
  double sigSDET() const { return sigTarg[4] - sigTarg[1] - sigTarg[2]; }

  // The target double diffractive excitation cross section.
  double sigDDE() const { return sigTarg[2]; }

  // The target non-diffractive (absorptive) cross section.
  double sigND() const { return sigTarg[1]; }

  // The Low-energy cross sections.
  double sigLow() const { return sigLExc() + sigLAnn() + sigLRes(); }

  // The Low-energy excitation cross section (code 157).
  double sigLExc() const { return sigTarg[8]; }

  // The Low-energy annihilation cross section (code 158).
  double sigLAnn() const { return sigTarg[9]; }

  // The Low-energy resonant cross section (code 159).
  double sigLRes() const { return sigTarg[10]; }

  // The target elastic b-slope parameter.
  double bSlope() const { return sigTarg[7]; }

  // Return the average non-diffractive impact parameter.
  double avNDB() const { return avNDb; }

  // Update internally stored cross sections.
  void updateSig(int idAIn, int idBIn, double eCMIn);

  // Calculate the Chi2 for the given cross section estimates.
  double Chi2(const SigEst & sigs, int npar) const;

  // Set beam kinematics.
  bool setKinematics(double eCMIn);

  // Set projectile particle.
  bool setIDA(int idA);

  // Use a genetic algorithm to fit the parameters.
  bool evolve(int nGenerations, double eCM, int idANow);

  // Get the number of free parameters for the model.
  int nParms() const { return parmSave.size(); }

  // Set the parameters of this model.
  void setParm(const vector<double>& parmIn) {
    for (size_t i = 0; i < parmSave.size(); ++i)
      parmSave[i] = parmIn[i];
  }

  // Get the current parameters of this model.
  vector<double> getParm() const { return parmSave; }

  // Get the minimum allowed parameter values for this model.
  virtual vector<double> minParm() const = 0;

  // Get the default parameter values for this model.
  virtual vector<double> defParm() const = 0;

  // Get the maximum allowed parameter values for this model.
  virtual vector<double> maxParm() const = 0;

  // Generate possible states for the nucleons in the projectile and
  // target nuclei.
  virtual void generateNucleonStates(Nucleus&, Nucleus&) {}

  // Take two nuclei and produce the corresponding subcollisions. The states
  // of the nucleons may be changed if fluctuations are allowed by the model.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) = 0;

  // Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const = 0;

private:

  // Generate parameters based on run settings and the evolutionary algorithm.
  bool genParms();

  // Save/load parameter configuration to/from settings/disk.
  bool saveParms(string fileName) const;
  bool loadParms(string fileName);

protected:

  // The nucleon-nucleon cross sections targets for this model
  // (tot, nd, dd, sdp, sdt, cd, el, bslope) and the required precision.
  vector<double> sigTarg, sigErr;
  vector< vector<double> > sigTargNN;

  // Saved parameters.
  vector<double> parmSave;

  // The parameters steering the fitting of internal parameters to
  // the different nucleon-nucleon cross sections.
  int NInt, NPop;
  double sigFuzz;
  double impactFudge;
  bool fitPrint;
  double eCMlow;

  // The estimated average impact parameter distance (in femtometer)
  // for absorptive collisions.
  double avNDb;
  double avNDolap;

  // Info from the controlling HeavyIons object.
  NucleusModel* projPtr          = {};
  NucleusModel* targPtr          = {};
  SigmaTotal* sigTotPtr          = {};
  SigmaCombined* sigCmbPtr       = {};
  Settings* settingsPtr          = {};
  Info* infoPtr                  = {};
  Rndm* rndmPtr                  = {};
  Logger* loggerPtr              = {};
  ParticleData * particleDataPtr = {};

  // For variable energies.
  int idASave, idBSave;
  bool doVarECM, doVarBeams;
  double eMin, eMax, eSave;
  int eCMPts;

  // The list of particles that have been fitted.
  vector<int> idAList;

  // A vector of interpolators for the current particle. Each entry
  // corresponds to one parameter, each interpolator is over the energy range.
  vector<LogInterpolator> *subCollParmsPtr;

  // Mapping id -> interpolator, one entry for each particle.
  map<int, vector<LogInterpolator>> subCollParmsMap;

  // Generation of elastic NN scatterings is turned on by default,
  // even if the cross section comes out wrong.
  int elasticMode;

  // Add in elastic NN scatterings to get the generated inelastic AA
  // cross section better described.
  double elasticFudge;

  // Helper class to cache cross sections at low energy
  struct SigmaCache {

    map<pair<int,int>, map<int, vector< vector<double> > > > cache;

    ParticleData * particleDataPtr = {};
    SigmaCombined * sigCmbPtr = {};
    double eCMStep = 0.1;

    void set(vector< vector<double> > & sigNN,
             int idAIn, int idBIn, double eCM) {
      auto & beamCache = cache[{idAIn, idBIn}];
      int iECMl = int(eCM/eCMStep);
      int iECMh = iECMl + 1;
      auto itl = beamCache.find(iECMl);
      if ( itl == beamCache.end() ) {
        fillCache(sigNN, idAIn, idBIn, iECMl);
        itl = beamCache.find(iECMl);
      }
      auto ith = beamCache.find(iECMh);
      if ( ith == beamCache.end() ) {
        fillCache(sigNN, idAIn, idBIn, iECMh);
        ith = beamCache.find(iECMh);
      }
      auto & sigl = itl->second;
      auto & sigh = ith->second;
      double fl = (eCM - iECMl*eCMStep)/eCMStep;
      double fh = 1.0 - fl;
      for ( size_t i = 0; i < sigl.size(); ++i )
        for ( size_t j = 0; j < sigl[i].size(); ++j )
          sigNN[i][j] =
            ( sigl[i][j] > 0.0? fl*sigl[i][j] + fh*sigh[i][j]: 0.0 );
    }

    void fillCache(vector< vector<double> > & sigNN,
                   int idAIn, int idBIn, int iECM) {
      double eCMIn = iECM*eCMStep;
      // Loop over both protons and neutrons.
      for ( int i = 0; i < 4; ++i ) {
        int idA = idAIn;
        int idB = idBIn;
        if ( i%2 == 1 ) {
          if ( abs(idA) != 2212 ) continue;
          idA = ( idA > 0? 2112: -2112 );
        }
        if ( i/2 == 1 ) {
          if ( abs(idB) != 2212 ) continue;
          idB = ( idB > 0? 2112: -2112 );
        }
        double mA = particleDataPtr->m0(idA);
        double mB = particleDataPtr->m0(idB);
        // Total.
        sigNN[i][0] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 0)*MB2FMSQ;
        // Non-diffractive.
        sigNN[i][1] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 1)*MB2FMSQ;
        // Doubly diffractive.
        sigNN[i][2] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 5)*MB2FMSQ;
        // Diffractive (and wounded) projectile.
        sigNN[i][3] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 4)*MB2FMSQ +
          sigNN[i][1] + sigNN[i][2];
        // Diffractive (and wounded) target.
        sigNN[i][4] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 3)*MB2FMSQ +
          sigNN[i][1] + sigNN[i][2];
        // Central diffractive.
        sigNN[i][5] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 6)*MB2FMSQ;
        // Elastic.
        sigNN[i][6] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 2)*MB2FMSQ;
        // b-slope not used for low energy.
        sigNN[i][7] = 0.0;
        // Low energy excitation.
        sigNN[i][8] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 7)*MB2FMSQ;
        // Low energy annihilation
        sigNN[i][9] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 8)*MB2FMSQ;
        // Low energy resonance
        sigNN[i][10] =
          sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 9)*MB2FMSQ;
        // Does not work for K0.
        if ( idA == 130 ) sigNN[i][10] = 0;
      }
      cache[{idAIn, idBIn}][iECM] = sigNN;
    }

  };

  SigmaCache lowEnergyCache;

};

//==========================================================================

// The most naive sub-collision model, assuming static nucleons and
// an absorptive cross section equal to the total inelastic. No
// fluctuations, meaning no diffraction.

class BlackSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  BlackSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~BlackSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> defParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Get cross sections used by this model.
  virtual SigEst getSig() const override {
    SigEst s;
    s.sig[0] = sigTot();
    s.sig[1] = sigND();
    s.sig[6] = s.sig[0] - s.sig[1];
    s.sig[7] = bSlope();
    return s;
  }

  // Take two nuclei and return the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A very simple sub-collision model, assuming static nucleons and
// just assuring that the individual nucleon-nucleon cross sections
// are preserved.

class NaiveSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  NaiveSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~NaiveSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> defParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Get cross sections used by this model.
  virtual SigEst getSig() const override {
    SigEst s;
    s.sig[0] = sigTot();
    s.sig[1] = sigND();
    s.sig[3] = sigSDEP();
    s.sig[4] = sigSDET();
    s.sig[2] = sigDDE();
    s.sig[6] = sigEl();
    s.sig[7] = bSlope();
    return s;
  }

  // Take two nuclei and return the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A base class for sub-collision models where each nucleon has a
// fluctuating "radius". The base model has two parameters, sigd and alpha,
// which are used for opacity calculations. Subclasses may have additional
// parameters to describe the radius distributions of that specific model.

class FluctuatingSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  FluctuatingSubCollisionModel(int nParmIn, int modein)
    : SubCollisionModel(nParmIn + 2), opacityMode(modein),
      sigd(parmSave[nParmIn]), alpha(parmSave[nParmIn + 1]) {}


  // Virtual destructor.
  virtual ~FluctuatingSubCollisionModel() override {}

  // Virtual init method.
  virtual bool init(int idAIn, int idBIn, double eCMIn) override;

  // Generate fluctuating radii for the nucleons in the projectile and
  // target nuclei.
  virtual void generateNucleonStates(Nucleus& proj, Nucleus& targ) override;

  // Take two nuclei and pick specific states for each nucleon,
  // then get the corresponding sub-collisions.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;
  SubCollisionSet getCollisionsNew(Nucleus& proj, Nucleus& targ);
  // Helper function.
  vector<double> getCollTypeProbs(const vector<double> & T) const;

  // Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const override;

protected:

  // Pick a radius for the nucleon, depending on the specific model.
  virtual double pickRadiusProj() const = 0;
  virtual double pickRadiusTarg() const = 0;

  // Optional mode for opacity.
  int opacityMode;

private:

  // Saturation scale of the nucleus.
  double& sigd;

  // Power of the saturation scale
  double& alpha;

  // The opacity of the collision at a given sigma.
  double opacity(double sig) const {
    sig /= sigd;
    if ( opacityMode )
      return pow(-expm1(-sig), alpha);
    return sig > numeric_limits<double>::epsilon() ?
      pow(-expm1(-1.0/sig), alpha) : 1.0;
  }

  // Get the minus the log of 1-T0 safely even when T0 close to 1.
  double log1mT0(double sig) const {
    double T0 = opacity(sig);
    if ( 1.0 - T0 > 1.0e6*numeric_limits<double>::epsilon() )
        return -log(1.0 - T0);
    if ( opacityMode == 1 )
      return sig/sigd -log(alpha);
    else
      return sigd/sig -log(alpha);
  }

  /// Calculate the overlap of two nucleons.
  double getOverlap(double b, double rt, double rp) const {
    double sig = M_PI*pow2(rp+rt);
    double T0  = opacity(sig);
    // *** TODO *** check that T0 is never 0.
    rt/=sqrt(2*T0);
    rp/=sqrt(2*T0);
    return (b>(rt+rp))? 0.0: 2.0*log1mT0(sig)/(2.0*T0 - pow2(T0));
  }

  // Return the elastic amplitude for a projectile and target state
  // and the impact parameter between the corresponding nucleons.
  double Tpt(const Nucleon::State & p,
             const Nucleon::State & t, double b) const {
    double sig = M_PI*pow2(p[0] + t[0]);
    double grey = opacity(sig);
    return sig/grey > b*b*2.0*M_PI? grey: 0.0;
  }

};

//==========================================================================

// A sub-collision model where each nucleon has a fluctuating
// "radius" according to a Strikman-inspired distribution.

class DoubleStrikmanSubCollisionModel : public FluctuatingSubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  DoubleStrikmanSubCollisionModel(int modeIn = 0)
    : FluctuatingSubCollisionModel(1, modeIn), k0(parmSave[0]) {}

  // Virtual destructor.
  virtual ~DoubleStrikmanSubCollisionModel() override {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return {  0.01,  1.0,  0.0  }; }
  vector<double> defParm() const override { return {  2.15, 17.24, 0.33 }; }
  vector<double> maxParm() const override {
    return { (opacityMode == 0? 60.00: 10.0), 60.0,
             (opacityMode == 0? 20.0: 2.0) }; }

protected:

  double pickRadiusProj() const override {
    double r =  rndmPtr->gamma(k0, r0());
    return (r < numeric_limits<double>::epsilon() ?
      numeric_limits<double>::epsilon() : r);
  }
  double pickRadiusTarg() const override {
    double r =  rndmPtr->gamma(k0, r0());
    return (r < numeric_limits<double>::epsilon() ?
      numeric_limits<double>::epsilon() : r);
  }

private:

  // The power in the Gamma distribution.
  double& k0;

  // Return the average radius deduced from other parameters and
  // the total cross section.
  double r0() const {
    return sqrt(sigTot() / (M_PI * (2.0 * k0 + 4.0 * k0 * k0)));
  }

};

//==========================================================================

// ImpactParameterGenerator is able to generate a specific impact
// parameter together with a weight such that a weighted average over
// any quantity X(b) corresponds to the infinite integral over d^2b
// X(b). This base class gives a Gaussian profile, d^2b exp(-b^2/2w^2).

class ImpactParameterGenerator {

public:

  // The default constructor takes a general width (in femtometers) as
  // argument.
  ImpactParameterGenerator() = default;

  // Virtual destructor.
  virtual ~ImpactParameterGenerator() {}

  // Virtual init method.
  virtual bool init();
  void initPtr(Info & infoIn, SubCollisionModel & collIn,
    NucleusModel & projIn, NucleusModel & targIn);

  // Return a new impact parameter and set the corresponding weight provided.
  virtual Vec4 generate(double & weight) const;

  // Return the scaling of the cross section used together with the
  // weight in generate() to obtain the cross section. This is by
  // default 1 unless forceUnitWeight is specified.
  virtual double xSecScale() const {
    return forceUnitWeight? M_PI*pow2(width()*cut): 2.0*M_PI *pow2(width());
  }

  // Set the width (in femtometers).
  void width(double widthIn) { widthSave = widthIn; }

  // Get the width.
  double width() const { return widthSave; }

  // Update width based on the associated subcollision and nucleus models.
  void updateWidth();

private:

  // The width of a distribution.
  double widthSave = 0.0;

  // The cut multiplied with widthSave to give the maximum allowed
  // impact parameter.
  double cut = 3.0;

  // Sample flat instead of with a Gaussian.
  bool forceUnitWeight = false;

protected:

  // Pointers from the controlling HeavyIons object.
  Info* infoPtr{};
  SubCollisionModel* collPtr{};
  NucleusModel* projPtr{};
  NucleusModel* targPtr{};
  Settings* settingsPtr{};
  Rndm* rndmPtr{};
  Logger* loggerPtr{};

};

//==========================================================================

// A sub-collision model where each nucleon fluctuates independently
// according to a log-normal distribution. Nucleons in the projectile and
// target may fluctuate according to different parameters, which is relevant
// e.g. for hadron-ion collisions with generic hadron species.

class LogNormalSubCollisionModel : public FluctuatingSubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  LogNormalSubCollisionModel(int modeIn = 0)
    : FluctuatingSubCollisionModel(4, modeIn),
    kProj(parmSave[0]), kTarg(parmSave[1]),
    rProj(parmSave[2]), rTarg(parmSave[3]) {}

  // Virtual destructor.
  virtual ~LogNormalSubCollisionModel() {}

  //virtual SigEst getSig() const override;

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override {
    return { 0.01, 0.01, 0.10, 0.10,  1.00, 0.00 }; }
  vector<double> defParm() const override {
    return { 1.00, 1.00, 0.54, 0.54, 17.24, 0.33 }; }
  vector<double> maxParm() const override {
    return { 2.00, 2.00, 4.00, 4.00, 20.00, 2.00 }; }

protected:

  double pickRadiusProj() const override { return pickRadius(kProj, rProj); }
  double pickRadiusTarg() const override { return pickRadius(kTarg, rTarg); }

private:

  // The standard deviation of each log-normal distribution.
  double& kProj;
  double& kTarg;

  // The mean radius of each nucleon.
  double& rProj;
  double& rTarg;

  double pickRadius(double k0, double r0) const {
    double logSig = log(M_PI * pow2(r0)) + k0 * rndmPtr->gauss();
    return sqrt(exp(logSig) / M_PI);
  }
};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HISubCollisionModel_H
