// FragmentationFlavZpT.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// StringFlav, StringZ and StringPT classes.

#include "Pythia8/FragmentationFlavZpT.h"

namespace Pythia8 {

//==========================================================================

// The StringFlav class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Offset for different meson multiplet id values.
const int StringFlav::mesonMultipletCode[6]
  = { 1, 3, 10003, 10001, 20003, 5};

// Clebsch-Gordan coefficients for baryon octet and decuplet are
// fixed once and for all, so only weighted sum needs to be edited.
// Order: ud0 + u, ud0 + s, uu1 + u, uu1 + d, ud1 + u, ud1 + s.
const double StringFlav::baryonCGOct[6]
  = { 0.75, 0.5, 0., 0.1667, 0.0833, 0.1667};
const double StringFlav::baryonCGDec[6]
  = { 0.,  0.,  1., 0.3333, 0.6667, 0.3333};

//--------------------------------------------------------------------------

// Initialize data members of the flavour generation.

void StringFlav::init() {

  // Set the fragmentation weights container.
  if (flag("VariationFrag:flav") || !infoPtr->weightContainerPtr
    ->weightsFragmentation.weightParms[WeightsFragmentation::Flav].empty())
    wgtsPtr = &infoPtr->weightContainerPtr->weightsFragmentation;

  // Basic parameters for generation of new flavour.
  probQQtoQ       = parm("StringFlav:probQQtoQ");
  probStoUD       = parm("StringFlav:probStoUD");
  probSQtoQQ      = parm("StringFlav:probSQtoQQ");
  probQQ1toQQ0    = parm("StringFlav:probQQ1toQQ0");

  // Spin parameters for combining two quarks to a diquark.
  vector<double> pQQ1tmp = settingsPtr->pvec("StringFlav:probQQ1toQQ0join");
  for (int i = 0; i < 4; ++i)
    probQQ1join[i] = 3. * pQQ1tmp[i] / (1. + 3. * pQQ1tmp[i]);

  // Parameters for normal meson production.
  for (int i = 0; i < 4; ++i) mesonRate[i][0] = 1.;
  mesonRate[0][1] = parm("StringFlav:mesonUDvector");
  mesonRate[1][1] = parm("StringFlav:mesonSvector");
  mesonRate[2][1] = parm("StringFlav:mesonCvector");
  mesonRate[3][1] = parm("StringFlav:mesonBvector");

  // Parameters for L=1 excited-meson production.
  mesonRate[0][2] = parm("StringFlav:mesonUDL1S0J1");
  mesonRate[1][2] = parm("StringFlav:mesonSL1S0J1");
  mesonRate[2][2] = parm("StringFlav:mesonCL1S0J1");
  mesonRate[3][2] = parm("StringFlav:mesonBL1S0J1");
  mesonRate[0][3] = parm("StringFlav:mesonUDL1S1J0");
  mesonRate[1][3] = parm("StringFlav:mesonSL1S1J0");
  mesonRate[2][3] = parm("StringFlav:mesonCL1S1J0");
  mesonRate[3][3] = parm("StringFlav:mesonBL1S1J0");
  mesonRate[0][4] = parm("StringFlav:mesonUDL1S1J1");
  mesonRate[1][4] = parm("StringFlav:mesonSL1S1J1");
  mesonRate[2][4] = parm("StringFlav:mesonCL1S1J1");
  mesonRate[3][4] = parm("StringFlav:mesonBL1S1J1");
  mesonRate[0][5] = parm("StringFlav:mesonUDL1S1J2");
  mesonRate[1][5] = parm("StringFlav:mesonSL1S1J2");
  mesonRate[2][5] = parm("StringFlav:mesonCL1S1J2");
  mesonRate[3][5] = parm("StringFlav:mesonBL1S1J2");

  // Store sum over multiplets for Monte Carlo generation.
  for (int i = 0; i < 4; ++i) mesonRateSum[i]
    = mesonRate[i][0] + mesonRate[i][1] + mesonRate[i][2]
    + mesonRate[i][3] + mesonRate[i][4] + mesonRate[i][5];

  // Parameters for uubar - ddbar - ssbar meson mixing.
  for (int spin = 0; spin < 6; ++spin) {
    double theta;
    if      (spin == 0) theta = parm("StringFlav:thetaPS");
    else if (spin == 1) theta = parm("StringFlav:thetaV");
    else if (spin == 2) theta = parm("StringFlav:thetaL1S0J1");
    else if (spin == 3) theta = parm("StringFlav:thetaL1S1J0");
    else if (spin == 4) theta = parm("StringFlav:thetaL1S1J1");
    else                theta = parm("StringFlav:thetaL1S1J2");
    double alpha = (spin == 0) ? 90. - (theta + 54.7) : theta + 54.7;
    alpha *= M_PI / 180.;
    // Fill in (flavour, spin)-dependent probability of producing
    // the lightest or the lightest two mesons of the nonet.
    mesonMix1[0][spin] = 0.5;
    mesonMix2[0][spin] = 0.5 * (1. + pow2(sin(alpha)));
    mesonMix1[1][spin] = 0.;
    mesonMix2[1][spin] = pow2(cos(alpha));
  }

  // Additional suppression of eta and etaPrime.
  etaSup      = parm("StringFlav:etaSup");
  etaPrimeSup = parm("StringFlav:etaPrimeSup");

  // Sum of baryon octet and decuplet weights.
  decupletSup = parm("StringFlav:decupletSup");
  for (int i = 0; i < 6; ++i) baryonCGSum[i]
    = baryonCGOct[i] + decupletSup * baryonCGDec[i];

  // Maximum SU(6) weight for ud0, ud1, uu1 types.
  baryonCGMax[0] = max( baryonCGSum[0], baryonCGSum[1]);
  baryonCGMax[1] = baryonCGMax[0];
  baryonCGMax[2] = max( baryonCGSum[2], baryonCGSum[3]);
  baryonCGMax[3] = baryonCGMax[2];
  baryonCGMax[4] = max( baryonCGSum[4], baryonCGSum[5]);
  baryonCGMax[5] = baryonCGMax[4];

  // Popcorn baryon parameters.
  popcornRate    = parm("StringFlav:popcornRate");
  popcornSpair   = parm("StringFlav:popcornSpair");
  popcornSmeson  = parm("StringFlav:popcornSmeson");

  // Suppression of leading (= first-rank) baryons.
  suppressLeadingB = flag("StringFlav:suppressLeadingB");
  lightLeadingBSup = parm("StringFlav:lightLeadingBSup");
  heavyLeadingBSup = parm("StringFlav:heavyLeadingBSup");

  // Enhanced-rate prefactor for MPIs and/or nearby string pieces.
  closePacking     = flag("ClosePacking:doClosePacking");
  enhanceStrange   = parm("ClosePacking:enhanceStrange");
  doEnhanceDiquark = flag("ClosePacking:doEnhanceDiquark");
  enhanceDiquark   = parm("ClosePacking:enhanceDiquark");
  exponentMPI      = parm("ClosePacking:expMPI");
  exponentNSP      = parm("ClosePacking:expNSP");

  // Save "vacuum" parameters for closepacking init() function.
  probStoUDSav    = probStoUD;
  probQQtoQSav    = probQQtoQ;
  probSQtoQQSav   = probSQtoQQ;
  probQQ1toQQ0Sav = probQQ1toQQ0;
  alphaQQSav      = (1. + 2. * probSQtoQQ * probStoUD + 9. * probQQ1toQQ0
    + 6. * probSQtoQQ * probQQ1toQQ0 * probStoUD
    + 3. * probQQ1toQQ0 * pow2(probSQtoQQ * probStoUD)) / (2. + probStoUD);

  // Calculate derived parameters.
  initDerived();

}

//--------------------------------------------------------------------------

// Initialise parameters when using close packing.

void StringFlav::init(double kappaModifier, double strangeJunc,
  double probQQmod) {

  double kappaRatio    = 1. + enhanceStrange * kappaModifier;
  double kappaInvRatio = 1. / pow(kappaRatio, 2*exponentNSP);

  // Altered probabilities with close packing.
  probStoUD    = pow(probStoUDSav, kappaInvRatio * (1 - strangeJunc));
  probSQtoQQ   = pow(probSQtoQQSav, kappaInvRatio);
  probQQ1toQQ0 = pow(probQQ1toQQ0Sav, kappaInvRatio);
  probQQtoQ    = probQQtoQSav;

  // If allowing effective kappa to enhance baryon production, do this.
  if (doEnhanceDiquark) {
    double alphaQQ = 1. + 2. * probSQtoQQ * probStoUD + 9. * probQQ1toQQ0
      + 6. * probSQtoQQ * probQQ1toQQ0 * probStoUD
      + 3. * probQQ1toQQ0 * pow2(probSQtoQQ * probStoUD);
    alphaQQ /= (2 + probStoUD);

    // Diquark scaling power controlled by enhanceDiquark.
    kappaRatio = 1. + enhanceDiquark * kappaModifier;
    kappaInvRatio   = 1. / pow(kappaRatio, 2*exponentNSP);
    probQQtoQ = alphaQQ * pow( (probQQtoQSav / alphaQQSav ), kappaInvRatio);
  }

  // Probability of a diquark being formed can scale with the probability
  // of a fluctuation on a string to not connect (and break) a nearby string.
  // for x probability of reconnection,
  // probability of diquark survival scales with 1/2 * [(1-x)^nG + (1-x)^nB]
  probQQtoQ = probQQmod * probQQtoQ;

  // Make sure probability is sensible.
  if (probQQtoQ > 1.) probQQtoQ = 1.;

  // Calculate derived parameters.
  initDerived();
}

//--------------------------------------------------------------------------

// Pick a new flavour (including diquarks) given an incoming one for
// Gaussian pTq^2 distribution.

FlavContainer StringFlav::pick(FlavContainer& flavOld, double, double,
  bool allowPop) {

  // Initial values for new flavour.
  FlavContainer flavNew;
  flavNew.rank = flavOld.rank + 1;

  // For original diquark assign popcorn quark and whether popcorn meson.
  int idOld = abs(flavOld.id);
  if (flavOld.rank == 0 && idOld > 1000 && allowPop) assignPopQ(flavOld);

  // Diquark exists, to be forced into baryon now.
  bool doOldBaryon    = (idOld > 1000 && flavOld.nPop == 0);
  // Diquark exists, but do meson now.
  bool doPopcornMeson = flavOld.nPop > 0;
  // Newly created diquark gives baryon now, antibaryon later.
  bool doNewBaryon    = false;

  // Choose whether to generate a new meson or a new baryon.
  if (!doOldBaryon && !doPopcornMeson && probQandQQ * rndmPtr->flat() > 1.) {
    doNewBaryon = true;
    if ((1. + popFrac) * rndmPtr->flat() > 1.) flavNew.nPop = 1;
  }

  // Optional suppression of first-rank baryon.
  if (flavOld.rank == 0 && doNewBaryon && suppressLeadingB) {
    double leadingBSup = (idOld < 4) ? lightLeadingBSup : heavyLeadingBSup;
    if (rndmPtr->flat() > leadingBSup) {
      doNewBaryon = false;
      flavNew.nPop = 0;
    }
  }

  // Single quark for new meson or for baryon where diquark already exists.
  if (!doPopcornMeson && !doNewBaryon) {
    flavNew.id = pickLightQ();
    if ( (flavOld.id > 0 && flavOld.id < 9) || flavOld.id < -1000 )
      flavNew.id = -flavNew.id;

    // Count breaks for variations, then done for simple-quark case.
    if (wgtsPtr != nullptr)
      wgtsPtr->flavStore(abs(flavNew.id), true, doOldBaryon);
    return flavNew;
  }

  // Case: 0 = q -> B B, 1 = q -> B M B, 2 = qq -> M B.
  int iCase = flavNew.nPop;
  if (flavOld.nPop == 1) iCase = 2;

  // Flavour of popcorn quark (= q shared between B and Bbar).
  if (doNewBaryon) {
    double sPopWT = dWT[iCase][0];
    if (iCase == 1) sPopWT *= scbBM[0] * popcornSpair;
    double rndmFlav = (2. + sPopWT) * rndmPtr->flat();
    flavNew.idPop = 1;
    if (rndmFlav > 1.) flavNew.idPop = 2;
    if (rndmFlav > 2.) flavNew.idPop = 3;
  } else flavNew.idPop = flavOld.idPop;

  // Flavour of vertex quark.
  double sVtxWT = dWT[iCase][1];
  if (flavNew.idPop >= 3) sVtxWT = dWT[iCase][2];
  if (flavNew.idPop > 3) sVtxWT *= 0.5 * (1. + 1./dWT[iCase][4]);
  double rndmFlav = (2. + sVtxWT) * rndmPtr->flat();
  flavNew.idVtx = 1;
  if (rndmFlav > 1.) flavNew.idVtx = 2;
  if (rndmFlav > 2.) flavNew.idVtx = 3;

  // Special case for light flavours, possibly identical.
  if (flavNew.idPop < 3 && flavNew.idVtx < 3) {
    flavNew.idVtx = flavNew.idPop;
    if (rndmPtr->flat() > dWT[iCase][3]) flavNew.idVtx = 3 - flavNew.idPop;
  }

  // Pick 2 * spin + 1.
  int spin = 3;
  if (flavNew.idVtx != flavNew.idPop) {
    double spinWT = dWT[iCase][6];
    if (flavNew.idVtx == 3) spinWT = dWT[iCase][5];
    if (flavNew.idPop >= 3) spinWT = dWT[iCase][4];
    if ((1. + spinWT) * rndmPtr->flat() < 1.) spin = 1;
  }

  // Form outgoing diquark. Count breaks for variations. Done.
  flavNew.id = 1000 * max(flavNew.idVtx, flavNew.idPop)
    + 100 * min(flavNew.idVtx, flavNew.idPop) + spin;
  if ( (flavOld.id < 0 && flavOld.id > -9) || flavOld.id > 1000 )
    flavNew.id = -flavNew.id;
  if (wgtsPtr != nullptr)
    wgtsPtr->flavStore(abs(flavNew.id), false, doOldBaryon);
  return flavNew;

}

//--------------------------------------------------------------------------

// Combine two flavours (including diquarks) to produce a hadron.
// The weighting of the combination may fail, giving output 0.

int StringFlav::combine(FlavContainer& flav1, FlavContainer& flav2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(flav1.id);
  int id2Abs = abs(flav2.id);
  int idMax  = max(id1Abs, id2Abs);
  int idMin  = min(id1Abs, id2Abs);

  // Construct a meson.
  if (idMax < 9 || idMin > 1000) {

    // Popcorn meson: use only vertex quarks. Fail if none.
    if (idMin > 1000) {
      id1Abs = flav1.idVtx;
      id2Abs = flav2.idVtx;
      idMax  = max(id1Abs, id2Abs);
      idMin  = min(id1Abs, id2Abs);
      if (idMin == 0) return 0;
    }

    // Pick spin state and preliminary code.
    int flav = (idMax < 3) ? 0 : idMax - 2;
    double rndmSpin = mesonRateSum[flav] * rndmPtr->flat();
    int spin = -1;
    do rndmSpin -= mesonRate[flav][++spin];
    while (rndmSpin > 0.);
    int idMeson = 100 * idMax + 10 * idMin + mesonMultipletCode[spin];

    // For nondiagonal mesons distinguish particle/antiparticle.
    if (idMax != idMin) {
      int sign = (idMax%2 == 0) ? 1 : -1;
      if ( (idMax == id1Abs && flav1.id < 0)
        || (idMax == id2Abs && flav2.id < 0) ) sign = -sign;
      idMeson *= sign;

    // For light diagonal mesons include uubar - ddbar - ssbar mixing.
    } else if (flav < 2) {
      double rMix = rndmPtr->flat();
      if      (rMix < mesonMix1[flav][spin]) idMeson = 110;
      else if (rMix < mesonMix2[flav][spin]) idMeson = 220;
      else                                   idMeson = 330;
      idMeson += mesonMultipletCode[spin];

      // Additional suppression of eta and eta' may give failure.
      if (idMeson == 221 && etaSup < rndmPtr->flat()) return 0;
      if (idMeson == 331 && etaPrimeSup < rndmPtr->flat()) return 0;
    }

    // Finished for mesons.
    return idMeson;
  }

  // SU(6) factors for baryon production may give failure.
  int idQQ1 = idMax / 1000;
  int idQQ2 = (idMax / 100) % 10;
  int spinQQ = idMax % 10;
  int spinFlav = spinQQ - 1;
  if (spinFlav == 2 && idQQ1 != idQQ2) spinFlav = 4;
  if (idMin != idQQ1 && idMin != idQQ2) spinFlav++;
  if (spinFlav < 0 || spinFlav > 5) return 0;
  if (baryonCGSum[spinFlav] < rndmPtr->flat() * baryonCGMax[spinFlav])
    return 0;

  // Order quarks to form baryon. Pick spin.
  int idOrd1 = max( idMin, max( idQQ1, idQQ2) );
  int idOrd3 = min( idMin, min( idQQ1, idQQ2) );
  int idOrd2 = idMin + idQQ1 + idQQ2 - idOrd1 - idOrd3;
  int spinBar = (baryonCGSum[spinFlav] * rndmPtr->flat()
    < baryonCGOct[spinFlav]) ? 2 : 4;

  // Distinguish Lambda- and Sigma-like.
  bool LambdaLike = false;
  if (spinBar == 2 && idOrd1 > idOrd2 && idOrd2 > idOrd3) {
    LambdaLike = (spinQQ == 1);
    if (idOrd1 != idMin && spinQQ == 1) LambdaLike = (rndmPtr->flat() < 0.25);
    else if (idOrd1 != idMin)           LambdaLike = (rndmPtr->flat() < 0.75);
  }

  // Form baryon code and return with sign.
  int idBaryon = (LambdaLike)
    ? 1000 * idOrd1 + 100 * idOrd3 + 10 * idOrd2 + spinBar
    : 1000 * idOrd1 + 100 * idOrd2 + 10 * idOrd3 + spinBar;
   return (flav1.id > 0) ? idBaryon : -idBaryon;

}

//--------------------------------------------------------------------------

// Combine three (di-)quark flavours into two hadrons.
// Note that at least one of the id's must be a diquark.

pair<int,int> StringFlav::combineDiquarkJunction(int id1, int id2, int id3) {

  // Order the junction ends in an increasing |id| sequence.
  if (abs(id1) > abs(id2)) swap(id1, id2);
  if (abs(id2) > abs(id3)) swap(id2, id3);
  if (abs(id1) > abs(id2)) swap(id1, id2);

  // If the first is a diquark then all are diquarks. Then split the first.
  // Combine its two quarks with the other diquarks into two baryons.
  int id1a = id1/1000;
  int id1b = (id1/100) % 10;
  int id2a = id2;
  int id2b = id3;

  // Otherwise the first is an antiquark. If the second is a diquark, also the
  // third is it. Then split the second. Let one of its quarks form a meson
  // with the first antiquark, and the other a baryon with the third diquark.
  if ( id1a == 0) {
    id1a = id2/1000;
    id1b = (id2/100) % 10;
    id2a = id1;
    id2b = id3;
  }

  // Finally, if the first two are antiquarks, then split the third diquark
  // and form two mesons.
  if ( id1a == 0) {
    id1a = id3/1000;
    id1b = (id3/100) % 10;
    id2a = id1;
    id2b = id2;
  }

  // If there was no diquark to split something is wrong.
  if (id1a == 0) return {0, 0};

  // Randomize the flavours of the split diquark and return the two hadrons.
  if (rndmPtr->flat() < 0.5) swap(id1a, id1b);
  return {combineId(id1a, id2a), combineId(id1b, id2b)};

}

//---------------------------------------------------------------------------

// Combine two flavours (including diquarks) to produce the lightest hadron
// allowed for that flavour content. No popcorn flavours.

int StringFlav::combineToLightest( int id1, int id2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int idMax  = max(id1Abs, id2Abs);
  int idMin  = min(id1Abs, id2Abs);
  int diqSgn = 0;

  // Quark-antiquark to meson.
  if (idMax < 9) {
    int idMeson = 100 * idMax + 10 * idMin + 1;

    // For nondiagonal mesons distinguish particle/antiparticle.
    if (idMax != idMin) {
      int sign = (idMax%2 == 0) ? 1 : -1;
      if (diqSgn != 0) sign *= diqSgn;
      else if ( (idMax == id1Abs && id1 < 0)
          || (idMax == id2Abs && id2 < 0) ) sign = -sign;
      idMeson *= sign;
    }

    // For light diagonal mesons pick pi0 or eta.
    else if (idMax <  3) idMeson = 111;
    else if (idMax == 3) idMeson = 221;

    // Finished for mesons.
    return idMeson;
  }

  // Quark-diquark to baryon
  int idQQ1  = idMax / 1000;
  int idQQ2  = (idMax / 100) % 10;
  int idOrd1 = max( idMin, max( idQQ1, idQQ2) );
  int idOrd3 = min( idMin, min( idQQ1, idQQ2) );
  int idOrd2 = idMin + idQQ1 + idQQ2 - idOrd1 - idOrd3;

  // Create baryon. Special cases with spin 3/2 and lambdalike.
  int idBaryon = 1000 * idOrd1 + 100 * idOrd2 + 10 * idOrd3 + 2;
  if (idOrd3 == idOrd1) idBaryon += 2;
  else if (idOrd2 != idOrd1 && idOrd3 != idOrd2)
    idBaryon = 1000 * idOrd1 + 100 * idOrd3 + 10 * idOrd2 + 2;

  // Finished for baryons.
  return (id1 > 0) ? idBaryon : -idBaryon;
}

//--------------------------------------------------------------------------

// Assign popcorn quark inside an original (= rank 0) diquark.

void StringFlav::assignPopQ(FlavContainer& flav) {

  // Safety check that intended to do something.
  int idAbs = abs(flav.id);
  if (flav.rank > 0 || idAbs < 1000) return;

  // Make choice of popcorn quark.
  int id1 = (idAbs/1000)%10;
  int id2 = (idAbs/100)%10;
  double pop2WT = 1.;
  if      (id1 == 3) pop2WT = scbBM[1];
  else if (id1 >  3) pop2WT = scbBM[2];
  if      (id2 == 3) pop2WT /= scbBM[1];
  else if (id2 >  3) pop2WT /= scbBM[2];
  // Agrees with Patrik code, but opposite to intention??
  flav.idPop = ((1. + pop2WT) * rndmPtr->flat() > 1.) ? id2 : id1;
  flav.idVtx = id1 + id2 - flav.idPop;

  // Also determine if to produce popcorn meson.
  flav.nPop = 0;
  double popWT = popS[0];
  if (id1 == 3) popWT = popS[1];
  if (id2 == 3) popWT = popS[2];
  if (idAbs%10 == 1) popWT *= sqrt(probQQ1toQQ0);
  if ((1. + popWT) * rndmPtr->flat() > 1.) flav.nPop = 1;

}

//--------------------------------------------------------------------------

// Combine two quarks to produce a diquark.
// Normally according to production composition, but nonvanishing idHad
// means diquark from known hadron content, so use SU(6) wave function.

int StringFlav::makeDiquark(int id1, int id2, int idHad) {

  // Initial values.
  int idMin = min( abs(id1), abs(id2));
  int idMax = max( abs(id1), abs(id2));
  int spin = 1;

  // Select spin of diquark formed from two valence quarks in proton.
  // (More hadron cases??)
  if (abs(idHad) == 2212 || abs(idHad) == 2112) {
    if (idMin == 1 && idMax == 2 && rndmPtr->flat() < 0.75) spin = 0;

  // Else select spin of diquark according to assumed spin-1 suppression.
  } else if (idMin != idMax) {
    if (rndmPtr->flat() > probQQ1join[min(idMax,5) - 2]) spin = 0;
  }

  // Combined diquark code.
  int idNewAbs = 1000 * idMax + 100 * idMin + 2 * spin + 1;
  return (id1 > 0) ? idNewAbs : -idNewAbs;

}

//--------------------------------------------------------------------------

// Initialise the derived parameters.

void StringFlav::initDerived() {

  // Parameters derived from init calls.
  probQandQQ      = 1. + probQQtoQ;
  probQandS       = 2. + probStoUD;
  probQandSinQQ   = 2. + probSQtoQQ * probStoUD;
  probQQ1corr     = 3. * probQQ1toQQ0;
  probQQ1corrInv  = 1. / probQQ1corr;
  probQQ1norm     = probQQ1corr / (1. + probQQ1corr);

  // Enumerate distinguishable diquark types (in diquark first is popcorn q).
  enum Diquark {ud0, ud1, uu1, us0, su0, us1, su1, ss1};

  // Maximum SU(6) weight by diquark type.
  barCGMax[ud0] = baryonCGMax[0];
  barCGMax[ud1] = baryonCGMax[4];
  barCGMax[uu1] = baryonCGMax[2];
  barCGMax[us0] = baryonCGMax[0];
  barCGMax[su0] = baryonCGMax[0];
  barCGMax[us1] = baryonCGMax[4];
  barCGMax[su1] = baryonCGMax[4];
  barCGMax[ss1] = baryonCGMax[2];

  // Diquark SU(6) survival = Sum_quark (quark tunnel weight) * SU(6).
  double dMB[8];
  dMB[ud0] = 2. * baryonCGSum[0] + probStoUD * baryonCGSum[1];
  dMB[ud1] = 2. * baryonCGSum[4] + probStoUD * baryonCGSum[5];
  dMB[uu1] = baryonCGSum[2] + (1. + probStoUD) * baryonCGSum[3];
  dMB[us0] = (1. + probStoUD) * baryonCGSum[0] + baryonCGSum[1];
  dMB[su0] = dMB[us0];
  dMB[us1] = (1. + probStoUD) * baryonCGSum[4] + baryonCGSum[5];
  dMB[su1] = dMB[us1];
  dMB[ss1] = probStoUD * baryonCGSum[2] + 2. * baryonCGSum[3];
  for (int i = 1; i < 8; ++i) dMB[i] = dMB[i] / dMB[0];

  // Tunneling factors for diquark production; only half a pair = sqrt.
  double probStoUDroot    = sqrt(probStoUD);
  double probSQtoQQroot   = sqrt(probSQtoQQ);
  double probQQ1toQQ0root = sqrt(probQQ1toQQ0);
  double qBB[8];
  qBB[ud1] = probQQ1toQQ0root;
  qBB[uu1] = probQQ1toQQ0root;
  qBB[us0] = probSQtoQQroot;
  qBB[su0] = probStoUDroot * probSQtoQQroot;
  qBB[us1] = probQQ1toQQ0root * qBB[us0];
  qBB[su1] = probQQ1toQQ0root * qBB[su0];
  qBB[ss1] = probStoUDroot * pow2(probSQtoQQroot) * probQQ1toQQ0root;

  // spin * (vertex factor) * (half-tunneling factor above).
  double qBM[8];
  qBM[ud1] = 3. * qBB[ud1];
  qBM[uu1] = 6. * qBB[uu1];
  qBM[us0] = probStoUD * qBB[us0];
  qBM[su0] = qBB[su0];
  qBM[us1] = probStoUD * 3. * qBB[us1];
  qBM[su1] = 3. * qBB[su1];
  qBM[ss1] = probStoUD * 6. * qBB[ss1];

  // Combine above two into total diquark weight for q -> B Bbar.
  for (int i = 1; i < 8; ++i) qBB[i] = qBB[i] * qBM[i];

  // Suppression from having strange popcorn meson.
  qBM[us0] *= popcornSmeson;
  qBM[us1] *= popcornSmeson;
  qBM[ss1] *= popcornSmeson;

  // Suppression for a heavy quark of a diquark to fit into a baryon
  // on the other side of popcorn meson: (0) s/u for q -> B M;
  // (1) s/u for rank 0 diquark su -> M B; (2) ditto for s -> c/b.
  double inf = numeric_limits<double>::infinity();
  double uNorm = 1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1];
  scbBM[0] = (2. * (qBM[su0] + qBM[su1]) + qBM[ss1]) / uNorm;
  scbBM[1] = qBM[us0] != 0 ? scbBM[0] * popcornSpair * qBM[su0]/qBM[us0] : inf;
  scbBM[2] = (1. + qBM[ud1]) * (2. + qBM[us0]) / uNorm;

  // Include maximum of Clebsch-Gordan coefficients.
  for (int i = 1; i < 8; ++i) dMB[i] *= qBM[i];
  for (int i = 1; i < 8; ++i) qBM[i] *= barCGMax[i] / barCGMax[0];
  for (int i = 1; i < 8; ++i) qBB[i] *= barCGMax[i] / barCGMax[0];

  // Popcorn fraction for normal diquark production.
  double qNorm = uNorm * popcornRate / 3.;
  double sNorm = scbBM[0] * popcornSpair;
  popFrac = qNorm * (1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1]
    + sNorm * (qBM[su0] + qBM[su1] + 0.5 * qBM[ss1])) / (1. +  qBB[ud1]
    + qBB[uu1] + 2. * (qBB[us0] + qBB[us1]) + 0.5 * qBB[ss1]);

  // Popcorn fraction for rank 0 diquarks, depending on number of s quarks.
  popS[0] = qNorm * qBM[ud1] / qBB[ud1];
  popS[1] = qBB[us1] != 0 && qBB[su1] != 0 ? qNorm * 0.5 * (qBM[us1] / qBB[us1]
    + sNorm * qBM[su1] / qBB[su1]) : inf;
  popS[2] = qBB[ss1] != 0 ? qNorm * sNorm * qBM[ss1] / qBB[ss1] : inf;

  // Recombine diquark weights to flavour and spin ratios. Second index:
  // 0 = s/u popcorn quark ratio.
  // 1, 2 = s/u ratio for vertex quark if popcorn quark is u/d or s.
  // 3 = q/q' vertex quark ratio if popcorn quark is light and = q.
  // 4, 5, 6 = (spin 1)/(spin 0) ratio for su, us and ud.

  // Case 0: q -> B B.
  dWT[0][0] = (2. * (qBB[su0] + qBB[su1]) + qBB[ss1])
    / (1. + qBB[ud1] + qBB[uu1] + qBB[us0] + qBB[us1]);
  dWT[0][1] = 2. * (qBB[us0] + qBB[us1]) / (1. + qBB[ud1] + qBB[uu1]);
  dWT[0][2] = qBB[su0] + qBB[su1] != 0 ? qBB[ss1]/(qBB[su0] + qBB[su1]) : inf;
  dWT[0][3] = qBB[uu1] / (1. + qBB[ud1] + qBB[uu1]);
  dWT[0][4] = qBB[su0] != 0 ? qBB[su1] / qBB[su0] : inf;
  dWT[0][5] = qBB[us0] != 0 ? qBB[us1] / qBB[us0] : inf;
  dWT[0][6] = qBB[ud1];

  // Case 1: q -> B M B.
  dWT[1][0] = (2. * (qBM[su0] + qBM[su1]) + qBM[ss1])
    / (1. + qBM[ud1] + qBM[uu1] + qBM[us0] + qBM[us1]);
  dWT[1][1] = 2. * (qBM[us0] + qBM[us1]) / (1. + qBM[ud1] + qBM[uu1]);
  dWT[1][2] = qBM[su0] + qBM[su1] != 0 ? qBM[ss1]/(qBM[su0] + qBM[su1]) : inf;
  dWT[1][3] = qBM[uu1] / (1. + qBM[ud1] + qBM[uu1]);
  dWT[1][4] = qBM[su0] != 0 ? qBM[su1] / qBM[su0] : inf;
  dWT[1][5] = qBM[us0] != 0 ? qBM[us1] / qBM[us0] : inf;
  dWT[1][6] = qBM[ud1];

  // Case 2: qq -> M B; diquark inside chain.
  dWT[2][0] = (2. * (dMB[su0] + dMB[su1]) + dMB[ss1])
    / (1. + dMB[ud1] + dMB[uu1] + dMB[us0] + dMB[us1]);
  dWT[2][1] = 2. * (dMB[us0] + dMB[us1]) / (1. + dMB[ud1] + dMB[uu1]);
  dWT[2][2] = dMB[su0] + dMB[su1] != 0 ? dMB[ss1]/(dMB[su0] + dMB[su1]) : inf;
  dWT[2][3] = dMB[uu1] / (1. + dMB[ud1] + dMB[uu1]);
  dWT[2][4] = dMB[su0] != 0 ? dMB[su1] / dMB[su0] : inf;
  dWT[2][5] = dMB[us0] != 0 ? dMB[us1] / dMB[us0] : inf;
  dWT[2][6] = dMB[ud1];

}

//==========================================================================

// Functions for the Lund symmetric FF: unnormalised, average, and RMSD.

//--------------------------------------------------------------------------

// The unnormalised Lund FF

double LundFFRaw(double z, double a, double b, double c, double mT2) {

  if (z <= 0. || z >= 1.) return 0.;
  return pow(1. - z, a) / pow(z, c) * exp(-b * mT2 / z);

}

//--------------------------------------------------------------------------

// Average, <z>, of Lund FF.
// Return values:
//   > 0. : <z>.
//    -1. : failed to compute normalisation.
//    -2. : failed to compute <z>.

double LundFFAvg(double a, double b, double mT2, double tol = 1.e-6) {

  // Checks whether the integration succeeded.
  bool check;

  // Fragmentation function dependent on only z (defined as a lambda function).
  function<double(double)> lundFF;

  // Get denominator (lundFF is function of only z, c = 1).
  lundFF = [=](double z) { return LundFFRaw(z, a, b, 1., mT2); };
  double denominator = 1.;
  check = integrateGauss(denominator, lundFF, 0., 1., tol);
  if (!check || denominator <= 0.) return -1.;

  // Get numerator (lundFF is function of only z, c = 0).
  lundFF = [=](double z) { return LundFFRaw(z, a, b, 0., mT2); };
  double numerator = 0.;
  check = integrateGauss(numerator, lundFF, 0., 1., tol);
  if (!check || numerator <= 0.) return -2.;

  // Done.
  return numerator / denominator;

}

//--------------------------------------------------------------------------

// RMSD(z) = sqrt(<z^2> - <z>^2) of Lund FF.
// Return values:
//   > 0. : rmsd
//    -1. : failed to compute normalisation.
//    -2. : failed to compute <z>.
//    -3. : failed to compute <z^2>.

double LundFFRms(double a, double b, double mT2, double tol = 1.e-6) {

  // Checks whether the integration succeeded.
  bool check;

  // Fragmentation function dependent on only z (defined as a lambda function).
  function<double(double)> lundFF;

  // Get denominator (lundFF is function of only z, c = 1).
  lundFF = [=](double z) { return LundFFRaw(z, a, b, 1., mT2); };
  double denominator = 1.;
  check = integrateGauss(denominator, lundFF, 0., 1., tol);
  if (!check || denominator <= 0.) return -1.;

  // Get first moment (lundFF is function of only z, c = 0).
  lundFF = [=](double z) { return LundFFRaw(z, a, b, 0., mT2); };
  double moment1 = 0.;
  check = integrateGauss(moment1, lundFF, 0., 1., tol);
  if (!check || moment1 <= 0.) return -2.;

  // Get second moment (lundFF is function of only z, c = -1).
  lundFF = [=](double z) { return LundFFRaw(z, a, b, -1., mT2); };
  double moment2 = 0.;
  check = integrateGauss(moment2, lundFF, 0., 1., tol);
  if (!check || moment2 <= 0.) return -3.;

  // Done.
  return sqrt(moment2 / denominator - pow2(moment1 / denominator));

}

//==========================================================================

// The StringZ class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// When a or c are close to special cases, default to these.
const double StringZ::CFROMUNITY = 0.01;
const double StringZ::AFROMZERO  = 0.02;
const double StringZ::AFROMC     = 0.01;

// Do not take exponent of too large or small number.
const double StringZ::EXPMAX     = 50.;

//--------------------------------------------------------------------------

// Initialize data members of the string z selection.
// Returns true if initialisation succeeded, false if failed.

bool StringZ::init() {

  // c and b quark masses.
  mc2           = pow2( particleDataPtr->m0(4));
  mb2           = pow2( particleDataPtr->m0(5));

  // Paramaters of Lund/Bowler symmetric fragmentation function.
  aLund         = parm("StringZ:aLund");
  bLund         = parm("StringZ:bLund");
  aExtraSQuark  = parm("StringZ:aExtraSQuark");
  aExtraDiquark = parm("StringZ:aExtraDiquark");
  rFactC        = parm("StringZ:rFactC");
  rFactB        = parm("StringZ:rFactB");
  rFactH        = parm("StringZ:rFactH");

  // Alternative parameterisation of Lund FF it terms of its average and
  // optionally rms and multiplicative factors for aDiquark and aStrange.
  if ( mode("StringZ:deriveLundPars") >= 1 ) {
    bool deriveA   = mode("StringZ:deriveLundPars") >= 2;
    bool deriveAQQ = mode("StringZ:deriveLundPars") >= 3;
    bool deriveAS  = mode("StringZ:deriveLundPars") >= 4;
    if (!deriveABLund( deriveA, deriveAQQ, deriveAS )) {
      loggerPtr->ABORT_MSG("derivation of Lund FF parameters failed");
      return false;
    }
  }

  // Use old or new behavior for aExtraSQuark and aExtraDiquark
  useOldAExtra  = flag("StringZ:useOldAExtra");

  // Flags and parameters of nonstandard Lund fragmentation functions.
  useNonStandC  = flag("StringZ:useNonstandardC");
  useNonStandB  = flag("StringZ:useNonstandardB");
  useNonStandH  = flag("StringZ:useNonstandardH");
  aNonC         = parm("StringZ:aNonstandardC");
  aNonB         = parm("StringZ:aNonstandardB");
  aNonH         = parm("StringZ:aNonstandardH");
  bNonC         = parm("StringZ:bNonstandardC");
  bNonB         = parm("StringZ:bNonstandardB");
  bNonH         = parm("StringZ:bNonstandardH");

  // Flags and parameters of Peterson/SLAC fragmentation function.
  usePetersonC  = flag("StringZ:usePetersonC");
  usePetersonB  = flag("StringZ:usePetersonB");
  usePetersonH  = flag("StringZ:usePetersonH");
  epsilonC      = parm("StringZ:epsilonC");
  epsilonB      = parm("StringZ:epsilonB");
  epsilonH      = parm("StringZ:epsilonH");

  // Parameters for joining procedure.
  stopM         = parm("StringFragmentation:stopMass");
  stopNF        = parm("StringFragmentation:stopNewFlav");
  stopS         = parm("StringFragmentation:stopSmear");

  // Parameters for reweighting.
  zHead         = parm("VariationFrag:zHead");
  posthoc       = flag("VariationFrag:z");

  // Set the fragmentation weights container.
  if (posthoc || !infoPtr->weightContainerPtr->
    weightsFragmentation.weightParms[WeightsFragmentation::Z].empty())
    wgtsPtr = &infoPtr->weightContainerPtr->weightsFragmentation;


  // All is well.
  return true;

}

//--------------------------------------------------------------------------

// Alternative parameterisation of the Lund function. Derive the bLund
// parameter given the average z for fixed a and mT2.

double StringZ::deriveBLund(double avgZ, double a, double mT2ref) {

  // Define lundFF as a function of only b, fixing a, and mT2 as parameters.
  auto lundFF = [=](double b) { return LundFFAvg(a, b, mT2ref); };

  // Solve for b and return.
  bNow = -1;
  bool check = brent(bNow, lundFF, avgZ, 0.0, 20.0, 1.e-7);
  return check ? bNow : -1;

}

//--------------------------------------------------------------------------

// Method to derive bLund and, optionally, aLund, aExtraDiquark,
// and aExtraSQuark, from:
//      avgZLund = <z(rho)>,
//      rmsZLund = sqrt( <z(rho)^2> - <z(rho)>^2),
//      facALundDiquark = (aLund + aExtraDiquark)/aLund,
//      facALundSQuark  = (aLund + aExtraStrange)/aLund,
// for reference (typical) values of the transverse mass mT.

bool StringZ::deriveABLund( bool deriveA, bool deriveAExtraDiquark,
  bool deriveAExtraSQuark ) {

  // Set up using reference mT2ref = mHad^2 + 2*sigmaPT^2 with mHad =
  // mRho, mK*, mp+ for light mesons, strange mesons, and baryons.
  double mRef        = particleDataPtr->m0(113);
  double mT2ref      = pow2(mRef) + 2.*pow2(parm("StringPT:sigma"));
  double mRefQQ      = particleDataPtr->m0(2212);
  double mT2refQQ    = pow2(mRefQQ)
    + 2.*pow2(parm("StringPT:sigma")*parm("StringPT:widthPreQQ0"));
  double mRefS       = particleDataPtr->m0(323);
  double mT2refS     = pow2(mRefS)
    + 2.*pow2(parm("StringPT:sigma")*parm("StringPT:widthPreStrange"));
  double avgZ        = parm("StringZ:avgZLund");
  double rmsZ        = parm("StringZ:rmsZLund");
  double facAQQ      = parm("StringZ:facALundDiquark");
  double facAS       = parm("StringZ:facALundSQuark");
  // Always use same starting point for derived parameters,
  // so that results are independent of previous settings/inits.
  double aNow        = (deriveA) ? 0.5 : parm("StringZ:aLund") ;
         bNow        = 1.0;
  double aExtraQQNow = parm("StringZ:aExtraDiquark");
  double aExtraSNow  = parm("StringZ:aExtraSQuark");

  // Debug output if requested.
  bool doReport = settingsPtr->mode("Print:verbosity") >= 3;
  if (doReport) {
    cout << scientific << setprecision(3) << setw(9)
         << "\n Deriving Lund FF parameter(s) with avgZ = " << avgZ;
    if (deriveA) cout << " rmsZ = " << rmsZ;
    else cout << " aLund = " << aNow;
    if (deriveAExtraDiquark) cout << " facADiquark = " << facAQQ;
    else cout << " aExtraDiquark = " << aExtraDiquark;
    if (deriveAExtraSQuark) cout << " facASQuark = " << facAS;
    else cout << " aExtraSQuark = " << aExtraSQuark;
    cout << endl;
  }

  // Simplest option: just derive bLund from requested avgZ.
  if ( !deriveA ) {
    if (doReport) {
      double avgZNow = LundFFAvg(aNow, bNow, mT2ref, 1.e-6);
      double rmsZNow = LundFFRms(aNow, bNow, mT2ref, 1.e-7);
      if (doReport) cout << fixed
                         << "   For aNow = " << aNow << " bNow = " << bNow
                         << ", got avgZNow = " << avgZNow
                         << " rmsZNow = " << rmsZNow << endl;
    }
    bNow = deriveBLund( avgZ, aNow, mT2ref);
    if (bNow < 0) {
      loggerPtr->ERROR_MSG("unable to converge on bLund");
      return false;
    }
  } else {
    // Derive both aLund and bLund from requested avgZ and rmsZ.
    bool accept  = false;
    double nLoop = 0;
    while (!accept) {
      if (++nLoop > 10000.) {
        loggerPtr->ERROR_MSG("maximum number of iterations exceeded");
        break;
      }
      const double TOLAVGZ = 1.e-5;
      const double TOLRMSZ = 1.e-5;
      double avgZNow   = LundFFAvg(aNow, bNow, mT2ref, 1.e-7);
      double rmsZNow   = LundFFRms(aNow, bNow, mT2ref, 1.e-7);
      if (doReport) cout << scientific << setprecision(3) << setw(9)
                         << "   For aNow = " << aNow << " bNow = " << bNow
                         << "  =>  avgZNow = " << avgZNow
                         << " rmsZNow = " << rmsZNow << endl;
      double deltaAvg = avgZNow - avgZ;
      double deltaRms = rmsZNow - rmsZ;

      // Take big steps in the beginning, then smaller ones.
      double step;
      if (nLoop < 500) step = 20.;
      else if (nLoop < 1000) step = 10.;
      else if (nLoop < 2000) step = 5.;
      else if (nLoop < 5000) step = 2.;
      else step = 1.;

      if ( abs(deltaRms) > TOLRMSZ ) {
        // First see if we can get the right RMS.
        aNow *= (1. + min(0.1,max(-0.1, step*deltaRms)));
        bNow = deriveBLund( avgZ, aNow, mT2ref);
        // Stop if we cannot possibly get a bigger width.
        if (aNow <= 0.001 && LundFFRms(aNow, bNow, mT2ref, 1.e-7)
          + 2*TOLRMSZ < rmsZ) {
          loggerPtr->ERROR_MSG("requested rmsZLund gave aLund < 0: "
            "forcing aLund = 0");
          aNow = 0.0;
          bNow = deriveBLund( avgZ, aNow, mT2ref);
          break;
        }
      }
      else if ( abs(deltaAvg) > TOLAVGZ ) {
        // Then get the right mean.
        aNow *= (1. + min(0.1,max(-0.1, step*deltaAvg)));
        bNow = deriveBLund( avgZ, aNow, mT2ref);
      }
      else accept = true;
    }
    // Check if method produced physical values.
    if (aNow < 0. || bNow < 0.) {
      loggerPtr->ERROR_MSG("unable to converge");
      return false;
    }
  }

  // Derive aExtraDiquark if requested.
  if (deriveAExtraDiquark) aExtraQQNow = (facAQQ - 1)*aNow;

  // Derive aExtraStrange if requested.
  if (deriveAExtraSQuark) aExtraSNow = (facAS - 1)*aNow;

  // Print out derived value(s).
  if ( !settingsPtr->flag("Print:quiet") ) {
    cout << "\n *-------  PYTHIA Derivation of Lund FF Parameters ----------"
      "------------------------------------------------------*" << endl;
    cout << fixed << setprecision(3) << " |\n | aLund = " << aNow
         << " & bLund = " << bNow << " GeV^-2 accepted";
    cout << "  (=> avgZ(rho) = " << setw(5)
         << LundFFAvg(aNow, bNow, mT2ref, 1.e-6)
         << " & rmsZ(rho) = " << setw(5)
         << LundFFRms(aNow, bNow, mT2ref, 1.e-6)
         << " for mTref = " << setw(5) << sqrt(mT2ref) << " GeV)" << endl;
    cout << fixed << setprecision(3) << " | aExtraSQuark  = " << aExtraSNow
         << "   (=> avgZ(K*) = " << setw(5)
         << LundFFAvg(aNow + aExtraSNow, bNow, mT2refS, 1.e-6)
         << " & rmsZ(K*) = " << setw(5)
         << LundFFRms(aNow + aExtraSNow, bNow, mT2refS, 1.e-6)
         << " for mTref = " << setw(5) << sqrt(mT2refS) << " GeV)" << endl;
    cout << fixed << setprecision(3) << " | aExtraDiquark = " << aExtraQQNow
         << "   (=> avgZ(p+) = " << setw(5)
         << LundFFAvg(aNow + aExtraQQNow, bNow, mT2refQQ, 1.e-6)
         << " & rmsZ(p+) = " << setw(5)
         << LundFFRms(aNow + aExtraQQNow, bNow, mT2refQQ, 1.e-6)
         << " for mTref = " << setw(5) << sqrt(mT2refQQ) << " GeV)" << endl;
    cout << " |\n *-------  End PYTHIA Derivation of Lund FF Parameters "
      "------------------------------------------------------------*" << endl;
  }

  // Set and check if derived bLund fell inside the nominal range.
  bool outOfRange = false;
  settingsPtr->parm("StringZ:bLund", bNow, false);
  if ( bNow != parm("StringZ:bLund") ) {
    // If outside nominal range, force so fits can see behaviour.
    outOfRange = true;
    settingsPtr->parm("StringZ:bLund", bNow, true);
  }

  // Set and check if derived aLund fell inside the nominal range.
  if ( deriveA ) {
    settingsPtr->parm("StringZ:aLund", aNow, false);
    if ( aNow != parm("StringZ:aLund") ) {
      // If outside nominal range, force so fits can see behaviour.
      outOfRange = true;
      settingsPtr->parm("StringZ:aLund", aNow, true);
    }
  }

  // Set and check if derived aExtraDiquark fell inside the nominal range.
  if ( deriveAExtraDiquark ) {
    settingsPtr->parm("StringZ:aExtraDiquark", aExtraQQNow, false);
    if ( aExtraQQNow != parm("StringZ:aExtraDiquark") ) {
      // If outside nominal range, force so fits can see behaviour.
      outOfRange = true;
      settingsPtr->parm("StringZ:aExtraDiquark", aExtraQQNow, true);
    }
  }

  // Set and check if derived aExtraDiquark fell inside the nominal range.
  if ( deriveAExtraSQuark ) {
    settingsPtr->parm("StringZ:aExtraSQuark", aExtraSNow, false);
    if ( aExtraSNow != parm("StringZ:aExtraSQuark") ) {
      // If outside nominal range, force so fits can see behaviour.
      outOfRange = true;
      settingsPtr->parm("StringZ:aExtraSQuark", aExtraSNow, true);
    }
  }

  // Issue warning if one or more parameters out of range.
  if (outOfRange) {
    loggerPtr->WARNING_MSG("one or more parameters out of range (forced)");
  }

  // No further calls needed since parameters updated in settings database.
  settingsPtr->mode("StringZ:deriveLundPars", 0);
  return true;

}

//--------------------------------------------------------------------------

// Initialize the flavour parameters.

void StringZ::initFlav(int idOld, int idNew) {

  // Find if old or new flavours correspond to diquarks.
  int idOldAbs = abs(idOld);
  int idNewAbs = abs(idNew);
  isOldSQuark = (idOldAbs == 3);
  isNewSQuark = (idNewAbs == 3);
  isOldDiquark = (idOldAbs > 1000 && idOldAbs < 10000);
  isNewDiquark = (idNewAbs > 1000 && idNewAbs < 10000);

  // Find heaviest quark in fragmenting parton/diquark.
  idFrag = idOldAbs;
  if (isOldDiquark) idFrag = max( idOldAbs / 1000, (idOldAbs / 100) % 10);

}

//--------------------------------------------------------------------------

// Initialize the shape parameters.

void StringZ::initShape(double mT2) {

  // Nonstandard a and b values implemented for heavy flavours.
  double aNow = aLund;
  bNow = bLund;
  if (idFrag == 4 && useNonStandC) {
    aNow = aNonC;
    bNow = bNonC;
  } else if (idFrag == 5 && useNonStandB) {
    aNow = aNonB;
    bNow = bNonB;
  } else if (idFrag >  5 && useNonStandH) {
    aNow = aNonH;
    bNow = bNonH;
  }

  // Shape parameters of Lund symmetric fragmentation function.
  aShape = aNow;
  // Old behavior used a_old instead of a_new in the
  // (1-z)^a factor for strange quarks and diquarks.
  // This is a bug but is kept for older tune compatibility.
  if (useOldAExtra) {
    if (isOldSQuark)  aShape += aExtraSQuark;
    if (isOldDiquark) aShape += aExtraDiquark;
  // This is the correct behavior that should by default be used.
  } else {
    if (isNewSQuark)  aShape += aExtraSQuark;
    if (isNewDiquark) aShape += aExtraDiquark;
  }
  bShape = bNow * mT2;
  cShape = 1.;
  if (isOldSQuark)  cShape -= aExtraSQuark;
  if (isNewSQuark)  cShape += aExtraSQuark;
  if (isOldDiquark) cShape -= aExtraDiquark;
  if (isNewDiquark) cShape += aExtraDiquark;
  if (idFrag == 4) cShape += rFactC * bNow * mc2;
  if (idFrag == 5) cShape += rFactB * bNow * mb2;
  if (idFrag >  5) cShape += rFactH * bNow * mT2;

}

//--------------------------------------------------------------------------

// Initialize the function sampling parameters.

void StringZ::initFunc(double a, double b, double c, double z, double zMax,
  double fPrel, double head) {

  bool aIsZero = (a < AFROMZERO);
  double fExp = b * (1. / zMax - 1. / z)+ c * log(zMax / z);
  if (!aIsZero) {
    if (z == 1) fExp = -numeric_limits<double>::infinity();
    else fExp += a * log( (1. - z) / (1. - zMax) );
  }
  fVal = exp( max( -EXPMAX, min( EXPMAX, fExp) ) ) ;
  fPrb = fVal / (fPrel * head);

}

//--------------------------------------------------------------------------

// Generate the fraction z that the next hadron will take,
// using either Lund/Bowler or, for heavy, Peterson/SLAC functions.
// Note: for a heavy new coloured particle we assume pT negligible.

double StringZ::zFrag(int idOld, int idNew, double mT2) {

  // Store the info needed for post-hoc reweighting.
  idNewNow = idNew;
  idOldNow = idOld;
  mT2Now   = mT2;

  // Use Peterson where explicitly requested for heavy flavours.
  initFlav(idOld, idNew);
  if (idFrag == 4 && usePetersonC) return zPeterson( epsilonC);
  if (idFrag == 5 && usePetersonB) return zPeterson( epsilonB);
  if (idFrag >  5 && usePetersonH) {
    double epsilon = epsilonH * mb2 / mT2;
    return zPeterson( epsilon);
  }

  // Determine the shape parameters and return the z.
  initShape(mT2);
  if (posthoc || (wgtsPtr != nullptr && !wgtsPtr->
      weightParms[WeightsFragmentation::Z].empty()))
    return zLund(aShape, bShape, cShape, zHead);
  else return zLund(aShape, bShape, cShape);

}

//--------------------------------------------------------------------------

// Determine the maximum for zLund.

double StringZ::zLundMax( double a, double b, double c) {

  // Normalization for Lund fragmentation function so that f <= 1.
  // Special cases for a = 0 and a = c.
  bool aIsZero = (a < AFROMZERO);
  bool aIsC = (abs(a - c) < AFROMC);

  // Determine position of maximum.
  double zMax;
  if (aIsZero) zMax = (c > b) ? b / c: 1.;
  else if (aIsC) zMax = b / (b + c);
  else { zMax = 0.5 * (b + c - sqrt( pow2(b - c) + 4. * a * b)) / (c - a);
    if (zMax > 0.9999 && b > 100.) zMax = min(zMax, 1. - a / b); }
  return zMax;
}

//--------------------------------------------------------------------------

// Generate a random z according to the Lund/Bowler symmetric
// fragmentation function f(z) = (1 -z)^a * exp(-b/z) / z^c.
// Normalized so that f(z_max) = 1  it can also be written as
// f(z) = exp( a * ln( (1 - z) / (1 - z_max) ) + b * (1/z_max - 1/z)
//           + c * ln(z_max/z) ).

// The arguments beginning with head are only needed for reweighting.

double StringZ::zLund(double a, double b, double c, double head) {

  // Special cases for c = 1, a = 0 and a = c.
  bool cIsUnity = (abs( c - 1.) < CFROMUNITY);
  bool aIsZero = (a < AFROMZERO);

  // Determine position of maximum.
  double zMax = zLundMax(a, b, c);

  // Subdivide z range if distribution very peaked near either endpoint.
  bool peakedNearZero = (zMax < 0.1);
  bool peakedNearUnity = (zMax > 0.85 && b > 1.);

  // Find integral of trial function everywhere bigger than f.
  // (Dummy start values.)
  double fIntLow = 1.;
  double fIntHigh = 1.;
  double fInt = 2.;
  double zDiv = 0.5;
  double zDivC = 0.5;
  // When z_max is small use that f(z)
  //   < 1     for z < z_div = 2.75 * z_max,
  //   < (z_div/z)^c for z > z_div (=> logarithm for c = 1, else power).
  if (peakedNearZero) {
    zDiv = 2.75 * zMax;
    fIntLow = zDiv;
    if (cIsUnity) fIntHigh = -zDiv * log(zDiv);
    else { zDivC = pow( zDiv, 1. - c);
           fIntHigh = zDiv * (1. - 1./zDivC) / (c - 1.);}
    fInt = fIntLow + fIntHigh;
  // When z_max large use that f(z)
  //   < exp( b * (z - z_div) ) for z < z_div with z_div messy expression,
  //   < 1   for z > z_div.
  // To simplify expressions the integral is extended to z =  -infinity.
  } else if (peakedNearUnity) {
    double rcb = sqrt(4. + pow2(c / b));
    zDiv = rcb - 1./zMax - (c / b) * log( zMax * 0.5 * (rcb + c / b) );
    if (!aIsZero) zDiv += (a/b) * log(1. - zMax);
    zDiv = min( zMax, max(0., zDiv));
    fIntLow = 1. / b;
    fIntHigh = 1. - zDiv;
    fInt = fIntLow + fIntHigh;
  }

  // Choice of z, preweighted for peaks at low or high z. (Dummy start values.)
  double z = 0.5;
  bool   accept = false;
  do {
    // Choice of z flat good enough for distribution peaked in the middle;
    // if not this z can be reused as a random number in general.
    z = rndmPtr->flat();
    double fPrel = 1.;
    // When z_max small use flat below z_div and 1/z^c above z_div.
    if (peakedNearZero) {
      if (fInt * rndmPtr->flat() < fIntLow) z = zDiv * z;
      else if (cIsUnity) {z = pow( zDiv, z); fPrel = zDiv / z;}
      else { z = pow( zDivC + (1. - zDivC) * z, 1. / (1. - c) );
             fPrel = pow( zDiv / z, c); }
    // When z_max large use exp( b * (z -z_div) ) below z_div
    // and flat above it.
    } else if (peakedNearUnity) {
      if (fInt * rndmPtr->flat() < fIntLow) {
        z = zDiv + log(z) / b;
        fPrel = exp( b * (z - zDiv) );
      } else z = zDiv + (1. - zDiv) * z;
    }

    // Evaluate actual f(z) (if in physical range) and correct.
    if (z > 0 && z < 1) {
      double fRnd = rndmPtr->flat();
      initFunc(a, b, c, z, zMax, fPrel, head);
      accept = fPrb > fRnd;

      // Loop over the variation parameters.
      if (wgtsPtr != nullptr) {
        if (posthoc)
          wgtsPtr->zStore(idOldNow, idNewNow, mT2Now, accept, z, fPrel);
        for (auto &parms : wgtsPtr->weightParms[WeightsFragmentation::Z]) {
          wgtsPtr->reweightValueByIndex(parms.second,
            wgtsPtr->zWeight(parms.first[0], parms.first[1], parms.first[2],
              parms.first[3], idOldNow, idNewNow,
              accept ? mT2Now : -mT2Now, z, fPrel));
        }
      }
    }
  } while (!accept);

  // Done.
  return z;

}

//--------------------------------------------------------------------------

// Generate a random z according to the Peterson/SLAC formula
// f(z) = 1 / ( z * (1 - 1/z - epsilon/(1-z))^2 )
//      = z * (1-z)^2 / ((1-z)^2 + epsilon * z)^2.

double StringZ::zPeterson( double epsilon) {

  double z;

  // For large epsilon pick z flat and reject,
  // knowing that 4 * epsilon * f(z) < 1 everywhere.
  if (epsilon > 0.01) {
    do {
      z = rndmPtr->flat();
      fVal = 4. * epsilon * z * pow2(1. - z)
        / pow2( pow2(1. - z) + epsilon * z);
    } while (fVal < rndmPtr->flat());
    return z;
  }

  // Else split range, using that 4 * epsilon * f(z)
  //   < 4 * epsilon / (1 - z)^2 for 0 < z < 1 - 2 * sqrt(epsilon)
  //   < 1                       for 1 - 2 * sqrt(epsilon) < z < 1
  double epsRoot = sqrt(epsilon);
  double epsComb = 0.5 / epsRoot - 1.;
  double fIntLow = 4. * epsilon * epsComb;
  double fInt = fIntLow + 2. * epsRoot;
  do {
    if (rndmPtr->flat() * fInt < fIntLow) {
      z = 1. - 1. / (1. + rndmPtr->flat() * epsComb);
      fVal = z * pow2( pow2(1. - z) / (pow2(1. - z) + epsilon * z) );
    } else {
      z = 1. - 2. * epsRoot * rndmPtr->flat();
      fVal = 4. * epsilon * z * pow2(1. - z)
        / pow2( pow2(1. - z) + epsilon * z);
    }
  } while (fVal < rndmPtr->flat());
  return z;

}

//==========================================================================

// The StringPT class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// To avoid division by zero one must have sigma > 0.
const double StringPT::SIGMAMIN     = 0.2;

//--------------------------------------------------------------------------

// Initialize data members of the string pT selection.

void StringPT::init() {

  // Parameters of the pT width and enhancements.
  double sigma     = parm("StringPT:sigma");
  sigmaQ           = sigma / sqrt(2.);
  enhancedFraction = parm("StringPT:enhancedFraction");
  enhancedWidth    = parm("StringPT:enhancedWidth");
  widthPreStrange  = parm("StringPT:widthPreStrange");
  widthPreQQ0      = parm("StringPT:widthPreQQ0");
  widthPreQQ1      = parm("StringPT:widthPreQQ1");
  useWidthPre      = (widthPreStrange != 1.) || (widthPreQQ0 != 1.)
    || (widthPreQQ1 != 1.);

  // Enhanced-width prefactor for MPIs and/or nearby string pieces.
  closePacking     = flag("ClosePacking:doClosePacking");
  enhancePT        = parm("ClosePacking:enhancePT");
  exponentMPI      = parm("ClosePacking:expMPI");
  exponentNSP      = parm("ClosePacking:expNSP");

  // Parameter for pT suppression in MiniStringFragmentation.
  sigma2Had        = 2. * pow2( max( SIGMAMIN, sigma) );

  // Parameters for reweighting.
  posthoc          = flag("VariationFrag:pT");

  // Set the fragmentation weights container.
  if (posthoc || !infoPtr->weightContainerPtr->
    weightsFragmentation.weightParms[WeightsFragmentation::PT].empty())
    wgtsPtr = &infoPtr->weightContainerPtr->weightsFragmentation;

}

//--------------------------------------------------------------------------

// Generate Gaussian pT such that <p_x^2> = <p_x^2> = sigma^2 = width^2/2,
// but with small fraction multiplied up to a broader spectrum.
// The missing first arguent is idIn, if a flavour depedence is desired.

pair<double, double> StringPT::pxy(int idIn, double kappaModifier) {

  // Normal (classical) width selection and factor for sigma variations.
  double sigma = sigmaQ;
  double mult  = 1.;
  if (rndmPtr->flat() < enhancedFraction) mult = enhancedWidth;

  // Optional prefactor for strange quarks and/or diquarks.
  if (useWidthPre && abs(idIn) >= 3) {
    if (abs(idIn) > 1000) {
      if (abs(idIn)%10 == 3) mult *= widthPreQQ1;
      else mult *= widthPreQQ0;
    }
    mult *= pow(widthPreStrange, particleDataPtr->nQuarksInCode(idIn, 3) );
  }

  // Enhanced-width prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    mult *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    double kappaRatio = 1. + enhancePT * kappaModifier;
    mult *= pow(max(1.0, kappaRatio), exponentNSP);
  }
  sigma *= mult;

  // Generate (p_x, p_y) pair.
  pair<double, double> gauss2 = rndmPtr->gauss2();

  // Calculate the weights from the variations.
  double pT2 = pow2(gauss2.first) + pow2(gauss2.second);
  if (wgtsPtr != nullptr) {
    if (posthoc)
      wgtsPtr->pTStore(pT2, mult);
    for (auto &parms : wgtsPtr->weightParms[WeightsFragmentation::PT]) {
      wgtsPtr->reweightValueByIndex(parms.second,
        wgtsPtr->pTWeight(parms.first[0], pT2, mult));
    }
  }

  // Return the result.
  return pair<double, double>(sigma * gauss2.first, sigma * gauss2.second);

}

//==========================================================================

} // end namespace Pythia8
