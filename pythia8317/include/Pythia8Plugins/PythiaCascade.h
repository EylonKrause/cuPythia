// PythiaCascade.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.
// Author: Torbjorn Sjostrand.

#ifndef Pythia8_PythiaCascade_H
#define Pythia8_PythiaCascade_H

#include "Pythia8/Pythia.h"

namespace Pythia8 {

//==========================================================================

// Wrapper class for Pythia usage in cosmic ray or detector cascades.
// Here it is assumed that the user wants to control the full
// evolution.  The complete event record is then kept somewhere
// separate from PYTHIA, and specific input on production and decay
// generates the next subevent.

// Kindly note non-negligible changes as of Pythia 8.316. This is related
// both to bug fixes in Angantyr, and to a better understanding of what
// it is supposed to do. The main difference is that elastic scatterings
// now are not simulated, it being assumed that they have no impact on
// the evolution of the cascade. As a consequence cross sections and
// event properties also change.

// Intended flow:

// - init sets up all generation, given a maximum energy. This may be
//   time-consuming, less so if MPI initialization data can be reused.

// - sigmaSetuphN prepares a possible collision for a given hadron by
//   calculating the hadron-nucleon inelastic cross section (if possible).
//   Moderately time-consuming.

// - sigmahA uses the hN cross section calculated by sigaSetuphN and
//   returns the hadron-ion inelastic cross section for a collision
//   with a specified nucleus. It can be called several times to cover
//   a mix of target nuclei, with minimal time usage.

// - nextColl performs the hadron-nucleus collision, as a sequences of
//   hadron-nucleon ones. Can be quite time-consuming.

// - nextDecay can be used anytime to decay a particle. Each
//   individual decay is rather fast, but there may be many of them.

// - After a nextColl three methods can give information about the
//   collision just generated:
//   * nCollisions() gives the number of hadron-nucleon inelastic
//     (sub)collisions, i.e. the number of wounded target nucleons;
//   * firstCollisionCode() gives the process number of the first
//     (sub)collision, as enumerated in the Process Selection section
//     of the html manual.
//   * firstCollisionMPI() gives the number of multiparton interactions
//     in the first (sub)collision.

// - stat() can be used at end of run to give a summary of error
//   messages.  Negligible time usage.

// - references to particleData() and rndm() can be used in the main
// - program.

// Note: when a hadron interacts with a medium particle, the latter is
// added to the event record. This program uses these additional
// status codes for target particles in the medium:
//  181: the first (or only) target nucleon in a collision.
//  182: nucleons from further subcollisions with the same target nucleus.

// Warning: the large boosts involved at the higher cosmic ray
// energies are not perfectly managed, which leads to non-negligible
// energy-momentum non-conservation. The worst (sub)events are caught
// and regenerated.

//--------------------------------------------------------------------------

class PythiaCascade {

public:

  // Default constructor; all setup is done in init().
  PythiaCascade() = default;

  //--------------------------------------------------------------------------

  // Initialize PythiaCascade for a given maximal incoming energy.

  // Hadrons below the kinetic energy threshold eKinMin will not be allowed
  // to interact with the medium. If reduced from the 0.3 GeV default then
  // some inelastic interactions are still allowed, while others fail
  // (but the run keeps going).

  // The enhanceSDtarget provides an enhancement factor for target-side
  // single diffraction, process codes 104 and 154, in subcollisions
  // after the first. This gives topologies more similar to Angantyr ones.
  // If 0 then no SD enhancement, if 1 only SD, with default in between.
  // A larger value means a smaller multiplicity, and vice versa, so this
  // is a quick way to obtain an envelope of hadronization uncertainties.

  // The initFile, by default "../share/Pythia8/setups/InitDefaultMPI.cmnd",
  // provides MPI initialization  data over a range of CM-frame energies,
  // where the upper edge sets the limit for allowed collisions.
  // Use main424.cc to rerun if you change any of the parameters that
  // affect the MPI rates, such as the pT0 parameter, its energy dependence,
  // the alpha_strong value, and the choice parton distributions.

  // Keep rapidDecays = false if you want to do every decay yourself,
  // else all particle types with a tau0 below smallTau0 will be
  // decayed immediately.  The tau0 units are mm/c, so default gives
  // c * tau0 = 1e-10 mm = 100 fm.  Note that time dilation effects are
  // not included, and they can be large, hence the low default value.
  // The Pythia default setup for decaying particles is obtained for
  // rapidDecays = true and smallTau0 = 1000.

  // With slowDecays true the mu+-, pi+-, K+- and K0L are allowed to decay,
  // as is expected in an atmospheric cascade but not for collider studies.
  // Note that rapidDecays and smallTau0 still determine how those decays
  // are handled, in the competition between interactions and decays.

  // Keep listFinalOnly = false if you want to get back the full event record,
  // else only the "final" particles of the collision or decay are returned.

  bool init( double eKinMinIn = 0.3, double enhanceSDtargetIn = 0.5,
    string initFile = "../share/Pythia8/setups/InitDefaultMPI.cmnd",
    bool rapidDecaysIn = false, double smallTau0In = 1e-10,
    bool slowDecays = true, bool listFinalOnlyIn = false) {

    // Store input for future usage.
    eKinMin         = eKinMinIn;
    enhanceSDtarget = enhanceSDtargetIn;
    rapidDecays     = rapidDecaysIn;
    smallTau0       = smallTau0In;
    listFinalOnly   = listFinalOnlyIn;

    // Proton and neutron masses.
    mp      = pythiaMain.particleData.m0(2212);
    mn      = pythiaMain.particleData.m0(2112);

    // Main Pythia object for managing the cascade evolution in a nucleus.
    // Can also do decays, but no hard processes.
    pythiaMain.readString("ProcessLevel:all = off");
    if (slowDecays) {
      pythiaMain.readString("13:mayDecay  = on");
      pythiaMain.readString("211:mayDecay = on");
      pythiaMain.readString("321:mayDecay = on");
      pythiaMain.readString("130:mayDecay = on");
    }
    pythiaMain.settings.flag("ParticleDecays:limitTau0", rapidDecays);
    pythiaMain.settings.parm("ParticleDecays:tau0Max", smallTau0);

    // Reduce statistics printout to relevant ones.
    pythiaMain.readString("Print:quiet = on");
    pythiaMain.readString("Stat:showProcessLevel = off");
    pythiaMain.readString("Stat:showPartonLevel = off");

    // Initialize. Return if failure.
    if (!pythiaMain.init()) return false;

    // Secondary Pythia object for individual collisions, or decays.
    // Reuse existing MPI initialization file. Failure if not found.
    if ( !pythiaColl.readString("include = " + initFile)) {
      cout << "\n Abort: failed to find or read MPI initialization file"
           << endl;
      return false;
    }

    // Variable incoming beam type and energy.
    pythiaColl.readString("Beams:allowVariableEnergy = on");
    pythiaColl.readString("Beams:allowIDAswitch = on");

    // Initialization eCM energy according to initFile.
    eCMMax = pythiaColl.settings.parm("Beams:eCMMaxMPI");
    pythiaColl.settings.parm("Beams:eCM", eCMMax);

    // Must use the soft and low-energy QCD processes, except elastic.
    pythiaColl.readString("SoftQCD:inelastic = on");
    pythiaColl.readString("LowEnergyQCD:inelastic = on");

    // Primary (single) decay to be done by pythiaColl, to circumvent
    // limitTau0.
    if (slowDecays) {
      pythiaColl.readString("13:mayDecay  = on");
      pythiaColl.readString("211:mayDecay = on");
      pythiaColl.readString("321:mayDecay = on");
      pythiaColl.readString("130:mayDecay = on");
    }

    // Secondary decays to be done by pythiaMain, respecting limitTau0.
    pythiaColl.readString("HadronLevel:Decay = off");

    // Reduce printout and relax energy-momentum conservation.
    // (Unusually large errors unfortunate consequence of large boosts.)
    pythiaColl.readString("Print:quiet = on");
    pythiaColl.readString("Check:epTolErr = 0.01");
    pythiaColl.readString("Check:epTolWarn = 0.0001");
    pythiaColl.readString("Check:mTolErr = 0.01");

    // Redure statistics printout to relevant ones.
    pythiaColl.readString("Stat:showProcessLevel = off");
    pythiaColl.readString("Stat:showPartonLevel = off");

    // Initialize and done.
    return pythiaColl.init();

  }

  //--------------------------------------------------------------------------

  // Calculate the average number of inelastic hadron-nucleon (hN) collisions
  // in a hadron-nucleus one, as a function of the inelastic hN cross section.
  // Interpolate if not in table, assuming <n> - 1 propto A^{2/3}.

  double nCollAvg(int A) {

    // Studied nuclei by A number, with offset and slope of <nColl>(sigma):
    // 1H, 2H, 4He, 9Be, 12C, 14N, 16O, 27Al, 40Ar, 56Fe, 63Cu, 84Kr,
    // 107Ag, 129Xe, 197Au, 208Pb.
    static const int nA = 16;
    static const int tabA[] = {
      1, 2, 4, 9, 12, 14, 16, 27, 40, 56, 63, 84, 107, 129, 197, 208};
    static const double tabOffset[] = {
      0.0000, 0.0510, 0.1164, 0.2036, 0.2328, 0.2520, 0.2624, 0.3190,
      0.3562, 0.3898, 0.3900, 0.3446, 0.3496, 0.3504, 0.3484, 0.3415 };
    static const double tabSlope[] = {
      0.0000,0.00187,0.00496, 0.0107, 0.0136, 0.0152, 0.0169, 0.0243,
      0.0314, 0.0385, 0.0415, 0.0506, 0.0581, 0.0644, 0.0806, 0.0830 };
    static const double tabSlopeLo[] = {
      0.0000,0.00361,0.00884, 0.0174, 0.0210, 0.0233, 0.0252, 0.0340,
      0.0418, 0.0496, 0.0524, 0.0600, 0.0668, 0.0727, 0.0873, 0.0893 };

    for (int i = 0; i < nA; ++i) {
      if (A == tabA[i]) {
        return min( 1. + tabSlopeLo[i] * sigmaNow,
          1. + tabOffset[i] + tabSlope[i] * sigmaNow);
      } else if (A < tabA[i]) {
        double nColl1 = min( tabSlopeLo[i - 1] * sigmaNow,
          tabOffset[i - 1] + tabSlope[i - 1] * sigmaNow);
        double nColl2 = min( tabSlopeLo[i] * sigmaNow,
          tabOffset[i] + tabSlope[i] * sigmaNow);
        double wt1 = double(tabA[i] - A) / double(tabA[i] - tabA[i - 1]);
        return 1. + wt1 * pow( A / tabA[i - 1], 2./3.) * nColl1
          + (1. - wt1) * pow( A / tabA[i], 2./3.) * nColl2;
      }
    }

    return numeric_limits<double>::quiet_NaN();
  }

  //--------------------------------------------------------------------------

  // Calculate the hadron-nucleon (proton) inelastic collision cross section.
  // Return false if not possible to find.

  bool sigmaSetuphN(int idNowIn, Vec4 pNowIn, double mNowIn) {

    // Cannot (or does not want to) handle low-energy hadrons.
    if (pNowIn.e() - mNowIn < eKinMin) return false;

    // Cannot handle hadrons above maximum energy set at initialization.
    eCMNow = (pNowIn + Vec4(0, 0, 0, mp)).mCalc();
    if (eCMNow > eCMMax) {
      logger.ERROR_MSG("too high energy");
      return false;
    }

    // Save incoming quantities for reuse in later methods.
    idNow = idNowIn;
    pNow  = pNowIn;
    mNow  = mNowIn;

    // Calculate hadron-nucleon inelastic cross section.
    // Check if cross section vanishes.
    sigmaNow = pythiaColl.getSigmaTotal(idNow, 2212, eCMNow, mNow, mp)
      - pythiaColl.getSigmaPartial(idNow, 2212, eCMNow, mNow, mp, 2);
    if (sigmaNow < 0.001) {
      if (eCMNow - mNow - mp > eKinMin)
        logger.WARNING_MSG("vanishing cross section");
      return false;
    }

    // Done.
    return true;

  }

  //--------------------------------------------------------------------------

  // Calculate the hadron-nucleus cross section for a given nucleon
  // number A, using hN cross section from sigmaSetuphN. Interpolate
  // where not (yet) available.

  double sigmahA(int A) {

    // Restrict to allowed range 1 <= A <= 208.
    if (A < 1 || A > 208) {
      logger.ERROR_MSG("A is outside of valid range (1 <= A <= 208)");
      return 0.;
    }

    // Correction factor for number of h-nucleon collisions per
    // h-nucleus one.
    double sigmahA = A * sigmaNow / nCollAvg(A);

    // Done.
    return sigmahA;

  }

  //--------------------------------------------------------------------------

  // Generate a collision, and return the event record.
  // Input (Z, A) of nucleus, and optionally collision vertex.

  Event& nextColl(int Znow, int Anow, Vec4 vNow = Vec4() ) {

    // References to the two event records. Clear main event record.
    Event& eventMain = pythiaMain.event;
    Event& eventColl = pythiaColl.event;
    eventMain.clear();
    codeFirstColl = 0;
    nMPIsave = 0;

    // Restrict to allowed range 1 <= A <= 208.
    if (Anow < 1 || Anow > 208) {
      logger.ERROR_MSG("A is outside of valid range (1 <= A <= 208)");
      return eventMain;
    }

    // Insert incoming particle in cleared main event record.
    eventMain.append(90,   -11, 0, 0, 1, 1, 0, 0, pNow, mNow);
    int iHad = eventMain.append(idNow, 12, 0, 0, 0, 0, 0, 0, pNow, mNow);
    eventMain[iHad].vProd(vNow);

    // Set up for collisions on a nucleus.
    int np      = Znow;
    int nn      = Anow - Znow;
    int sizeOld = 0;
    int sizeNew = 0;
    Vec4 dirNow = pNow / pNow.pAbs();
    Rndm& rndm  = pythiaMain.rndm;

    // Drop rate of geometric series. (Deuterium is special case.)
    double probMore = (Anow == 2) ? nCollAvg(Anow) - 1.
                    : 1. - 1. / nCollAvg(Anow);
    nCollAcc    = 0;

    // Loop over varying number of hit nucleons in target nucleus.
    for (int iColl = 1; iColl <= Anow; ++iColl) {
      if (iColl > 1 && rndm.flat() > probMore) break;

      // Pick incoming projectile: trivial for first subcollision, else ...
      int iProj    = iHad;

      // ... find highest-pLongitudinal particle from latest subcollision.
      if (iColl > 1) {
        iProj = 0;
        double pMax = 0.;
        for (int i = sizeOld; i < sizeNew; ++i)
        if ( eventMain[i].isFinal() && eventMain[i].isHadron()) {
          double pp = dot3(dirNow, eventMain[i].p());
          if (pp > pMax) {
            iProj = i;
            pMax  = pp;
          }
        }

        // No further subcollision if no particle with enough energy.
        if ( iProj == 0
          || eventMain[iProj].e() - eventMain[iProj].m() < eKinMin) break;
      }

      // Pick one p or n from target.
      int idProj = eventMain[iProj].id();
      bool doProton = rndm.flat() < (np / double(np + nn));
      if (doProton) np -= 1;
      else          nn -= 1;
      int idNuc = (doProton) ? 2212 : 2112;

      // Current subcollision four-vectors and CM energy.
      Vec4 pProj = eventMain[iProj].p();
      double mTarg = (doProton) ? mp : mn;
      Vec4 pTarg( 0., 0., 0., mTarg);
      double eCMPT = (pProj + pTarg).mCalc();

      // Reject if process is only possible because projectile is off-shell.
      // (Current kinematics handling assumes that "beam particles" are
      // on the mass shell, but this could be changed eventually.)
      mProjMax = max(eventMain[iProj].m(), pythiaMain.particleData.m0(idProj));
      if (sqrt(pProj.pAbs2() + pow2(mProjMax)) - mProjMax < eKinMin) break;

      // Do a projectile-nucleon subcollision. Return empty event if failure.
      // Optionally enhance secondary single diffractive on target side.
      // for non-first collision, i.e. make more of codes 104 or 154.
      pythiaColl.setBeamIDs(idProj, idNuc);
      pythiaColl.setKinematics(eCMPT);
      int codeNow = (iColl > 1 && rndm.flat() < enhanceSDtarget) ? 4 : 0;
      if (!pythiaColl.next(codeNow)) {
        eventMain.clear();
        return eventMain;
      }

      // Statistics.
      if (iColl == 1) codeFirstColl = pythiaColl.info.code();
      if (iColl == 1) nMPIsave = pythiaColl.info.nMPI();

      // Boost back collision to lab frame.
      RotBstMatrix MtoLab;
      MtoLab.fromCMframe( pProj, pTarg);
      eventColl.rotbst( MtoLab);

      // Insert target nucleon. Mothers are (0,iProj) to mark who it
      // interacted with. Always use proton mass for simplicity.
      int statusNuc = (iColl == 1) ? -181 : -182;
      int iNuc = eventMain.append( idNuc, statusNuc, 0, iProj, 0, 0, 0, 0,
        0., 0., 0., mp, mp);
      eventMain[iNuc].vProdAdd(vNow);

      // Update full energy of the event with the target mass.
      eventMain[0].e( eventMain[0].e() + mTarg);
      eventMain[0].m( eventMain[0].p().mCalc() );

      // Insert secondary produced particles (but skip intermediate partons)
      // into main event record and shift to correct production vertex.
      sizeOld = eventMain.size();
      for (int iSub = 3; iSub < eventColl.size(); ++iSub) {
        if (!eventColl[iSub].isFinal()) continue;
        int iNew = eventMain.append(eventColl[iSub]);
        eventMain[iNew].mothers(iNuc, iProj);
        eventMain[iNew].vProdAdd(vNow);
      }
      sizeNew = eventMain.size();

      // Update daughters of colliding hadrons and other history.
      eventMain[iProj].daughters(sizeOld, sizeNew - 1);
      eventMain[iNuc].daughters(sizeOld, sizeNew - 1);
      eventMain[iProj].statusNeg();
      eventMain[iProj].tau(0.);

      // End of loop over interactions in a nucleus.
      ++nCollAcc;
    }

    // Optionally do decays of short-lived particles.
    if (rapidDecays) pythiaMain.moreDecays();

    // Optionally compress event record.
    if (listFinalOnly) compress();

    // Return generated collision.
    return eventMain;

  }

  //--------------------------------------------------------------------------

  // Generate a particle decay, and return the event record.
  // You can allow sequential decays, if they occur rapidly enough.

  Event& nextDecay(int idNowIn, Vec4 pNowIn, double mNowIn,
    Vec4 vNow = Vec4() ) {

    // Save incoming quantities. (Not needed, but by analogy with
    // collisions.)
    idNow = idNowIn;
    pNow  = pNowIn;
    mNow  = mNowIn;

    // References to the event records. Clear them.
    Event& eventMain = pythiaMain.event;
    Event& eventColl = pythiaColl.event;
    eventMain.clear();
    eventColl.clear();

    // Insert incoming particle in cleared collision event record.
    eventColl.append(90,   -11, 0, 0, 1, 1, 0, 0, pNow, mNow);
    int iHad = eventColl.append(idNow, 12, 0, 0, 0, 0, 0, 0, pNow, mNow);
    eventColl[iHad].vProd(vNow);

    // Decay incoming particle. Return empty event if fail. Copy event record.
    if (!pythiaColl.moreDecays(iHad)) return eventMain;
    eventMain = eventColl;

    // Optionally do secondary decays of short-lived particles.
    if (rapidDecays) pythiaMain.moreDecays();

    // Optionally compress event record.
    if (listFinalOnly) compress();

    // Return generated collision.
    return eventMain;

  }

  //--------------------------------------------------------------------------

  // Compress the event record by removing initial and intermediate
  // particles.  Keep line 0, since the += operator for event records
  // only adds from 1 on.

  void compress() {

    // Reference to the main event record. Original and new size.
    Event& eventMain = pythiaMain.event;
    int sizeOld = eventMain.size();
    int sizeNew = 1;

    // Loop through all particles and move up the final ones to the top.
    // Remove history information. Update event record size.
    for (int i = 1; i < sizeOld; ++i) if (eventMain[i].isFinal()) {
      eventMain[sizeNew] = eventMain[i];
      eventMain[sizeNew].mothers( 0, 0);
      eventMain[sizeNew].daughters( 0, 0);
      ++sizeNew;
    }

    // Shrink event record to new size.
    eventMain.popBack( sizeOld - sizeNew);

  }

  //--------------------------------------------------------------------------

  // Summary of aborts, errors and warnings.

  void stat() {
    pythiaMain.stat();
    pythiaColl.stat();
    logger.errorStatistics();
  }

  //--------------------------------------------------------------------------

  // Provide number of subcollisions, hardest subprocess, and whether elastic.
  // The latter should always be false after the recent code changes.

  int nCollisions() {return nCollAcc;}

  int firstCollisionCode() {return codeFirstColl;}

  int firstCollisionMPI() {return nMPIsave;}

  //--------------------------------------------------------------------------

  // Possibility to access particle data and random numbers from
  // pythiaMain.

  ParticleData& particleData() {return pythiaMain.particleData;}

  Rndm& rndm() {return pythiaMain.rndm;}

//--------------------------------------------------------------------------

private:

  // The Pythia instances used for decays and for collisions. Could
  // be made public, but cleaner to allow only limited access as
  // above.
  Pythia pythiaMain, pythiaColl;

  // Logger instance for errors in this class.
  Logger logger;

  // Save quantities.
  bool   rapidDecays, listFinalOnly;
  int    idNow, nCollAcc, codeFirstColl, nMPIsave;
  double eKinMin, enhanceSDtarget, smallTau0, mp, mn, eCMMax, mNow, mProjMax,
         eCMNow, sigmaNow;
  Vec4   pNow;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_PythiaCascade_H
