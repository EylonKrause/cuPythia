// ThermalFragmentation.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// ThermalFlav and ThermalPT classes.

#include "Pythia8/ThermalFragmentation.h"

namespace Pythia8 {

//==========================================================================

// The ThermalStringFlav class.

//--------------------------------------------------------------------------

// Initialize data members of the flavour generation.

void ThermalStringFlav::init() {

  // Temperature parameters for thermal model.
  temperature      = parm("StringPT:temperature");
  tempPreFactor    = parm("StringPT:tempPreFactor");

  // Hadron multiplets in thermal model.
  mesonNonetL1     = flag("StringFlav:mesonNonetL1");
  nNewQuark        = mode("StringFlav:nQuark");

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
    // Fill in rates for multiplication.
    mesMixRate1[0][spin] = mesonMix1[0][spin];
    mesMixRate2[0][spin] = mesonMix2[0][spin] - mesonMix1[0][spin];
    mesMixRate3[0][spin] = 1.0 - mesMixRate1[0][spin] - mesMixRate2[0][spin];
    mesMixRate1[1][spin] = mesonMix1[1][spin];
    mesMixRate2[1][spin] = mesonMix2[1][spin] - mesonMix1[1][spin];
    mesMixRate3[1][spin] = 1.0 - mesMixRate1[1][spin] - mesMixRate2[1][spin];
  }

  // Fill list of possible hadrons that are allowed to be produced.
  // Also include a list of "emergency" hadrons that are needed to get
  // rid of all possible endpoint (di)quarks.
  vector<int> hadIDsProd, hadIDsHvyC, hadIDsHvyB;

  // Baryon octet and decuplet.
  int baryonLight[18] = { 2112, 2212, 3112, 3122, 3212, 3222, 3312, 3322,
                          1114, 2114, 2214, 2224, 3114, 3214, 3224, 3314,
                          3324, 3334 };
  int baryonHvyC[22]  = { 4112, 4122, 4132, 4212, 4222, 4232, 4312, 4322,
                          4332, 4412, 4422, 4432,
                          4114, 4214, 4224, 4314, 4324, 4334, 4414, 4424,
                          4434, 4444 };
  int baryonHvyB[35]  = { 5112, 5122, 5132, 5142, 5212, 5222, 5232, 5242,
                          5312, 5322, 5332, 5342, 5412, 5422, 5432, 5442,
                          5512, 5522, 5532, 5542,
                          5114, 5214, 5224, 5314, 5324, 5334, 5414, 5424,
                          5434, 5444, 5514, 5524, 5534, 5544, 5554 };
  for (int i = 0; i < 18; i++) hadIDsProd.push_back( baryonLight[i] );
  // Check how many heavy baryons to include.
  if (nNewQuark > 4) {
    for (int i = 0; i < 35; i++) hadIDsProd.push_back( baryonHvyB[i] );
  } else {
    // Only include lightest combinations.
    int bBar[9] = { 5112, 5122, 5132, 5212, 5222, 5232, 5312, 5322, 5332 };
    for (int i = 0; i < 9; i++) {
      hadIDsHvyB.push_back(  bBar[i] );
      hadIDsHvyB.push_back( -bBar[i] );
    }
  }
  if (nNewQuark > 3) {
    for (int i = 0; i < 22; i++) hadIDsProd.push_back( baryonHvyC[i] );
  } else {
    // Only include lightest combinations.
    int cBar[9] = { 4112, 4122, 4132, 4212, 4222, 4232, 4312, 4322, 4332 };
    for (int i = 0; i < 9; i++) {
      hadIDsHvyC.push_back(  cBar[i] );
      hadIDsHvyC.push_back( -cBar[i] );
    }
  }
  // Antibaryons.
  int sizeNow = int(hadIDsProd.size());
  for (int i = 0; i < sizeNow; i++) hadIDsProd.push_back( -hadIDsProd[i] );

  // Mesons nonets. Take pseudoscalar PDG codes as basis.
  int mesonPSLight[9] = { 311, 321, 211, -311, -321, -211, 111, 221, 331 };
  int mesonPSHvyC[7]  = { 411, 421, 431, -411, -421, -431, 441 };
  int mesonPSHvyB[9]  = { 511, 521, 531, 541, -511, -521, -531, -541, 551 };
  vector<int> mesonPS;
  for (int i = 0; i < 9; i++) mesonPS.push_back( mesonPSLight[i] );
  // Check how many heavy mesons to include. If not included in ordinary
  // production, fill minimal list with "emergency" hadrons
  if (nNewQuark > 4) {
    for (int i = 0; i < 9; i++) mesonPS.push_back( mesonPSHvyB[i] );
  } else {
    // Include all possible combinations, only pseudoscalar as they
    // are the lightest ones.
    int bMes[10] = { 511, 521, 531, 541, -511, -521, -531, -541, 551 };
    for (int i = 0; i < 10; i++) hadIDsHvyB.push_back( bMes[i] );
  }
  if (nNewQuark > 3) {
    for (int i = 0; i < 7; i++) mesonPS.push_back( mesonPSHvyC[i] );
  } else {
    // Include all possible combinations, only pseudoscalar as they
    // are the lightest ones.
    int cMes[8] = { 411, 421, 431, -411, -421, -431, 441 };
    for (int i = 0; i < 8; i++) hadIDsHvyC.push_back( cMes[i] );
  }
  int nMeson = int(mesonPS.size());
  // Pseudoscalar nonet J=0, S=0, L=0.
  for (int i = 0; i < nMeson; i++)
    hadIDsProd.push_back( mesonPS[i] );
  // Vector nonet J=1, S=1, L=0.
  for (int i = 0; i < nMeson; i++)
    hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 2 : -2) );
  // Include L=1 nonets?
  if (mesonNonetL1) {
    // Pseudovector nonet J=1, S=0, L=1.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 10002 : -10002) );
    // Scalar nonet J=0, S=1, L=1.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 10000 : -10000) );
    // Pseudovector nonet J=1, S=1, L=1.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 20002 : -20002) );
    // Tensor nonet J=2, S=1, L=1.
    for (int i = 0; i < nMeson; i++)
      hadIDsProd.push_back( mesonPS[i] + (mesonPS[i] > 0 ? 4 : -4) );
  }

  // Fill list of all hadrons ids (ordinary and "emergency").
  vector<int> hadIDsAll;
  for (int i = 0; i < int(hadIDsProd.size()); i++)
    hadIDsAll.push_back( hadIDsProd[i] );
  for (int i = 0; i < int(hadIDsHvyC.size()); i++)
    hadIDsAll.push_back( hadIDsHvyC[i] );
  for (int i = 0; i < int(hadIDsHvyB.size()); i++)
    hadIDsAll.push_back( hadIDsHvyB[i] );

  // Fill map with IDs of hadron constituents for all hadrons.
  for (int i = 0; i < int(hadIDsAll.size()); i++) {
    int id    = hadIDsAll[i];
    int idAbs = abs(id);
    vector< pair<int,int> > quarkCombis;
    // Baryon can be split into q + qq in several different ways.
    if (particleDataPtr->isBaryon(id)) {
      bool isOctet   = ( (idAbs % 10) == 2 );
      int  q3        = (idAbs/10)   % 10;
      int  q2        = (idAbs/100)  % 10;
      int  q1        = (idAbs/1000) % 10;
      bool threeFlav = q1 != q2 && q1 != q3 && q2 != q3;
      // Baryon octet J=1/2.
      if (isOctet) {
        if (threeFlav) {
          // Add (q2+q3)_0/1 + q1.
          // if (q2 < q3) (q2+q3)_0 and if (q2 > q3) (q2+q3)_1.
          int j = (q2 < q3) ? 1 : 3;
          int qn[2]  = { min( q3, q2), max( q3, q2) };
          addQuarkDiquark(quarkCombis, q1,
            1000 * qn[1] + 100 * qn[0] + j, id);
          // Add other combinations. Can be both, J=0 or J=1.
          for (j = 1; j < 4; j += 2) {
            // (q1+q3)j + q2
            addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + j, id);
            // (q1+q2)j + q3
            addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + j, id);
          }
        } else {
          // Quarks with the same flavour form J=1,
          // all other combinations can be both, J=0 or J=1.
          for (int j = 1; j < 4; j += 2) {
            // (q1+q2)1 + q3
            if ( j == 3 || q1 != q2 )
              addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + j, id);
            // (q1+q3)1 + q2
            if ( j == 3 || q1 != q3 )
              addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + j, id);
            // (q2+q3)1 + q1
            if ( j == 3 || q2 != q3 )
              addQuarkDiquark(quarkCombis, q1, 1000 * q2 + 100 * q3 + j, id);
          }
        }
      }
      // Baryon decuplet J=3/2.
      else {
        // All quark pairs form diquarks with J=1.
        // (q1+q2)1 + q3
        addQuarkDiquark(quarkCombis, q3, 1000 * q1 + 100 * q2 + 3, id);
        // (q1+q3)1 + q2
        addQuarkDiquark(quarkCombis, q2, 1000 * q1 + 100 * q3 + 3, id);
        // (q2+q3)1 + q1
        addQuarkDiquark(quarkCombis, q1, 1000 * q2 + 100 * q3 + 3, id);
      }
    // Mesons usually have a trivial subdivision into quark + antiquark.
    // Mixing of diagonal mesons is taken into account later.
    } else {
      int q1        = (idAbs/100) % 10;
      bool uptype1  = (q1 % 2 == 0);
      int q2        = (idAbs/10)  % 10;
      bool uptype2  = (q2 % 2 == 0);
      int quark     = q1;
      int antiQuark = q2;
      // id > 0: downtype+uptype: up = quark, down = antiquark (default)
      // id > 0: same type -> larger id decides
      if ( uptype2 && !uptype1 ) swap( quark, antiQuark);
      if ( (q1 > q2 && !uptype1 && !uptype2)
        || (q2 > q1 &&  uptype2 &&   uptype1) ) swap( quark, antiQuark);
      if (id < 0) swap( quark, antiQuark);
      quarkCombis.push_back( make_pair( quark, -antiQuark) );
    }
    hadronConstIDs[id] = quarkCombis;
  }

  // Copy into smaller versions (one for ordinary production, two for
  // "emergency")
  map< int, vector< pair<int,int> > > hadConstIDsC, hadConstIDsB,
                                      hadConstIDsProd;
  for (int i=0; i<int(hadIDsAll.size()); i++) {
    int id = hadIDsAll[i];
    if (find(hadIDsProd.begin(), hadIDsProd.end(), id) != hadIDsProd.end())
      hadConstIDsProd[id] = hadronConstIDs[id];
    if (find(hadIDsHvyC.begin(), hadIDsHvyC.end(), id) != hadIDsHvyC.end())
      hadConstIDsC[id] = hadronConstIDs[id];
    if (find(hadIDsHvyB.begin(), hadIDsHvyB.end(), id) != hadIDsHvyB.end())
      hadConstIDsB[id] = hadronConstIDs[id];
  }
  map< int, map< int, vector< pair<int,int> > > > hadConstIDsHvy;
  hadConstIDsHvy[4] = hadConstIDsC;
  hadConstIDsHvy[5] = hadConstIDsB;

  // List with all possible initial (di)quarks we could get.
  int inIDs[26]    = { 1, 2, 3, 4, 5, 1103, 2203, 3303, 2101, 2103, 3101,
                       3103, 3201, 3203, 4101, 4103, 4201, 4203, 4301,
                       4303, 5101, 5103, 5201, 5203, 5301, 5303 };
  int inIDsHvyC[2] = { 4403, -4403 };
  int inIDsHvyB[6] = { 5503, -5503, 5401, -5401, 5403, -5403 };
  vector<int> incomingIDs;
  for (int i = 0; i < 26; i++) {
    incomingIDs.push_back( inIDs[i]);
    incomingIDs.push_back(-inIDs[i]);
  }
  // If we include heavy quark hadrons we include the following diquarks in
  // addition.
  if (nNewQuark > 3) {
    for (int i = 0; i < 2; i++) incomingIDs.push_back(inIDsHvyC[i]);
    if (nNewQuark > 4) {
      for (int i = 0; i < 6; i++) incomingIDs.push_back( inIDsHvyB[i]);
    }
  }
  int nIncome = int(incomingIDs.size());

  // Loop over list with all possible initial (di)quarks.
  // Fill map possibleHadrons with
  // key = initial (di)quark id, value = list of possible hadron ids
  //                                     + nr in hadronConstIDs.
  for (int iIDin = 0; iIDin < nIncome; iIDin++) {
    int idIn    = incomingIDs[iIDin];
    int idInAbs = abs(idIn);
    map< int, vector< pair<int,int> > > hadConstIDsNow = hadConstIDsProd;
    // For heavy quarks add "emergency" list, if needed.
    for (int iHvy = nNewQuark+1; iHvy <= 5; iHvy++) {
      if (particleDataPtr->nQuarksInCode(idInAbs, iHvy) > 0)
        for (map< int, vector< pair<int,int> > >::iterator
             it = hadConstIDsHvy[iHvy].begin();
             it != hadConstIDsHvy[iHvy].end(); ++it)
          hadConstIDsNow[it->first] = it->second;
    }
    // Fill list: first parameter of pair is hadron ID, second is nr of
    // hadron constituents in the list.
    vector< pair<int,int> > possibleHadronIDs;
    // Loop through list with hadrons and their (di)quark content,
    // check if possible to produce given the choice of initial (di)quark.
    for (map< int, vector< pair<int,int> > >::iterator
         it = hadConstIDsNow.begin(); it != hadConstIDsNow.end(); ++it) {
      vector< pair<int,int> > constituentIDs = it->second;
      int nConst   = int(constituentIDs.size());
      int hadronID = it->first;
      // Loop over constituent IDs.
      for (int iConst = 0; iConst < nConst; iConst++) {
        int ID1 = constituentIDs[iConst].first;
        int ID2 = constituentIDs[iConst].second;
        if ( (ID1 == idIn) || (ID2 == idIn) ) {
          possibleHadronIDs.push_back( make_pair(hadronID,iConst) );
          // To include uubar-ddbar-ssbar mixing include all diagonal mesons.
          if ( (idInAbs < 4) && (ID1 == -ID2) ) {
            if (idInAbs == 1) {
              possibleHadronIDs.push_back( make_pair(hadronID+110,iConst) );
              possibleHadronIDs.push_back( make_pair(hadronID+220,iConst) );
            } else if (idInAbs == 2) {
              possibleHadronIDs.push_back( make_pair(hadronID-110,iConst) );
              possibleHadronIDs.push_back( make_pair(hadronID+110,iConst) );
            } else if (idInAbs == 3) {
              possibleHadronIDs.push_back( make_pair(hadronID-220,iConst) );
              possibleHadronIDs.push_back( make_pair(hadronID-110,iConst) );
            }
          }
        }
      }
    }
    if (int(possibleHadronIDs.size()) < 1)
      loggerPtr->ERROR_MSG("no possible hadrons found");
    possibleHadrons[idIn] = possibleHadronIDs;
  }

  // Calculate baryon octet and decuplet weighting factors
  // based on Clebsch-Gordan coefficients and spin counting.
  // Parameters: qDi1 qDi2 q3 spin.
  // Zero for flavour=0 and same flavour diquarks with J=0.
  for (int q1 = 0; q1 < 6; q1++) {
    for (int q2 = 0; q2 < 6; q2++) {
      baryonOctWeight[q1][q1][q2][0] = 0.0; // qq0 + r
      baryonDecWeight[q1][q1][q2][0] = 0.0; // qq0 + r
      for (int spin = 0; spin < 1; spin++) {
        baryonOctWeight[ 0][q1][q2][spin] = 0.0;
        baryonOctWeight[q1][ 0][q2][spin] = 0.0;
        baryonOctWeight[q1][q2][ 0][spin] = 0.0;
        baryonDecWeight[ 0][q1][q2][spin] = 0.0;
        baryonDecWeight[q1][ 0][q2][spin] = 0.0;
        baryonDecWeight[q1][q2][ 0][spin] = 0.0;
      }
    }
  }
  // Clebsch-Gordon for the rest.
  for (int q1 = 1; q1 < 6; q1++) {
    baryonOctWeight[q1][q1][q1][1] = 0.0; // qq1 + q
    baryonDecWeight[q1][q1][q1][1] = 1.0;
    for (int q2 = 1; q2 < 6; q2++) if (q1!=q2) {
      baryonOctWeight[q1][q1][q2][1] = 0.1667; // qq1 + r
      baryonDecWeight[q1][q1][q2][1] = 0.3333;
      baryonOctWeight[q1][q2][q1][0] = 0.75;   // qr0 + q
      baryonDecWeight[q1][q2][q1][0] = 0.0;
      baryonOctWeight[q2][q1][q1][0] = 0.75;   // rq0 + q
      baryonDecWeight[q2][q1][q1][0] = 0.0;
      baryonOctWeight[q1][q2][q1][1] = 0.0833; // qr1 + q
      baryonDecWeight[q1][q2][q1][1] = 0.6667;
      baryonOctWeight[q2][q1][q1][1] = 0.0833; // rq1 + q
      baryonDecWeight[q2][q1][q1][1] = 0.6667;
      for (int q3 = 0; q3 < 6; q3++) if ((q1 != q3) && (q2 != q3)) {
        baryonOctWeight[q1][q2][q3][0] = 0.5;    // qr0 + s
        baryonDecWeight[q1][q2][q3][0] = 0.0;
        baryonOctWeight[q1][q2][q3][1] = 0.1667; // qr1 + s
        baryonDecWeight[q1][q2][q3][1] = 0.3333;
      }
    }
  }
  // Spin 1 diquarks get extra factor of 3. And all factors
  // get relative baryon-to-meson ratio.
  double BtoMratio = parm("StringFlav:BtoMratio");
  for (int q1 = 0; q1 < 6; q1++) {
    for (int q2 = 0; q2 < 6; q2++) {
      for (int q3 = 0; q3 < 6; q3++) {
        for (int spin = 0; spin < 2; spin++) {
          baryonOctWeight[q1][q2][q3][spin] *= BtoMratio;
          baryonDecWeight[q1][q2][q3][spin] *= BtoMratio;
          if (spin == 1) {
            baryonOctWeight[q1][q2][q3][1] *= 3.0;
            baryonDecWeight[q1][q2][q3][1] *= 3.0;
          }
        }
      }
    }
  }

  // Go through the list of possible hadrons and calculate the prefactor
  // that will multiply the rate.
  double strSup = parm("StringFlav:StrangeSuppression");
  for (int iIDin = 0; iIDin < nIncome; iIDin++) {
    int idIn      = incomingIDs[iIDin];
    int idInAbs   = abs(idIn);
    vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idIn];
    vector<double> prefactors;
    for (int iHad = 0; iHad < int(possibleHadronsNow.size()); iHad++) {
      double prefacNow = 1.0;
      // Get hadron and constituents.
      int hadronID    = possibleHadronsNow[iHad].first;
      int hadronIDabs = abs(hadronID);
      int iConst      = possibleHadronsNow[iHad].second;
      int ID1         = hadronConstIDs[hadronID][iConst].first;
      int ID2         = hadronConstIDs[hadronID][iConst].second;
      // Extra suppression factor for s/c/b quarks.
      double nHeavy   = 0.0;
      for (int i = 3; i <= 5; i++) {
        nHeavy += particleDataPtr->nQuarksInCode( ID1, i);
        nHeavy += particleDataPtr->nQuarksInCode( ID2, i);
      }
      prefacNow      *= pow(strSup, nHeavy);
      if (particleDataPtr->isMeson(hadronID)) {
        // Extra factor according to last digit for spin counting.
        prefacNow *= (abs(hadronID) % 10);
        // Include correct uubar-ddbar-ssbar mixing factor;
        if ( (idInAbs < 4) && (ID1 == -ID2) ) {
          int flav = ( (idInAbs < 3) ? 0 : 1 );
          // Get spin used as counter for the different multiplets
          int spin = getMesonSpinCounter(hadronID);
          double mesonMix[3] = { mesMixRate1[flav][spin],
                                 mesMixRate2[flav][spin],
                                 mesMixRate3[flav][spin] };
          prefacNow *= mesonMix[abs(ID1)-1];
        }
      } else {
        // Check if baryon is octet or decuplet.
        bool isOct = ((hadronIDabs % 10) == 2);
        // Make sure ID2 is diquark.
        if (abs(ID2) < abs(ID1)) swap(ID1,ID2);
        // Extract quark flavours and spin from diquark.
        int Q1 = ( (abs(ID2)/1000) % 10 );
        int Q2 = ( (abs(ID2)/100)  % 10 );
        if (Q1 > 5 || Q2 > 5) {
          loggerPtr->ERROR_MSG("invalid quark content flavours for diquark");
          continue;
        }
        int diqSpin = ( ((abs(ID2) % 10) == 1) ? 0 : 1 );
        // Single quark.
        int Q3      = abs(ID1);
        // Find Clebsch-Gordan: q1 in DQ | q2 in DQ | q3 | S of DQ
        if (isOct) prefacNow *= baryonOctWeight[Q1][Q2][Q3][diqSpin];
        else       prefacNow *= baryonDecWeight[Q1][Q2][Q3][diqSpin];
        // Special cases for Lamda (312) and Sigma (321) or the like.
        if ( isOct && (Q1!=Q2) && (Q1!=Q3) && (Q2!=Q3) ) {
          // Extract the two lightest quarks from hadron.
          int Qhad1   = ( (hadronIDabs/10)  % 10 );
          int Qhad2   = ( (hadronIDabs/100) % 10 );
          int QhadMin = min(Qhad1,Qhad2);
          int QhadMax = max(Qhad1,Qhad2);
          // Extract the two quarks from the diquark.
          int QdiqMin = min(Q1,Q2);
          int QdiqMax = max(Q1,Q2);
          // Don't do anything if (12) or (21) is diquark.
          if ( !((QdiqMin == QhadMin) && (QdiqMax == QhadMax)) ) {
            // Sigma (321)
            if (Qhad2 > Qhad1) prefacNow *= ( (diqSpin == 0) ? 0.75 : 0.25 );
            // Lamda (312)
            else               prefacNow *= ( (diqSpin == 0) ? 0.25 : 0.27 );
          }
        }
      }
      // Save prefactor.
      prefactors.push_back(prefacNow);
    }
    possibleRatePrefacs[idIn] = prefactors;
  }

  // Now the same again for joining the last two (di)quarks into hadron.
  for (int iIDin1 = 0; iIDin1 < nIncome; iIDin1++) {
    int idIn1     = incomingIDs[iIDin1];
    int idIn1Abs  = abs(idIn1);
    // Loop over possible partners, start with next quark.
    for (int iIDin2 = iIDin1+1; iIDin2 < nIncome; iIDin2++) {
      int idIn2      = incomingIDs[iIDin2];
      int idIn2Abs   = abs(idIn2);
      int idInNow[2] = { min(idIn1,idIn2), max(idIn1,idIn2) };
      pair<int,int> inPair = pair<int,int>(idInNow[0], idInNow[1]);
      // Skip all combinations with two diquarks.
      if ( (idIn1Abs > 1000) && (idIn2Abs > 1000) ) continue;
      // Skip all combinations with two quarks or two antiquarks.
      if ( ( ((idIn1 > 0) && (idIn2 > 0)) || ((idIn1 < 0) && (idIn2 < 0)) )
           && (idIn1Abs < 10) && (idIn2Abs < 10) ) continue;
      // Skip all combinations with quark-antidiquark and
      // antiquark-diquark. (1 = diquark, 2 = quark not possible).
      if ( ((idIn2 >  1000) && (idIn1Abs < 10) && (idIn1 < 0)) ||
           ((idIn2 < -1000) && (idIn1Abs < 10) && (idIn1 > 0)) ) continue;
      // If we are not including heavy quarks skip combinations
      // of heavy quark - diquark with heavy quark.
      if ((idIn1Abs < 10) && (idIn2Abs > 1000)) {
        vector< pair<int,int> > hvyCombs;
        if (nNewQuark < 5) {
          hvyCombs.push_back(make_pair(4,4));
          if (nNewQuark < 4) {
            hvyCombs.push_back(make_pair(5,4));
            hvyCombs.push_back(make_pair(4,5));
            hvyCombs.push_back(make_pair(5,5));
          }
        }
        bool skip = false;
        for (int iComb = 0; iComb < int(hvyCombs.size()); iComb++) {
          int idNow[2] = { hvyCombs[iComb].first, hvyCombs[iComb].second };
          if ( (particleDataPtr->nQuarksInCode(idIn2Abs,idNow[0]) > 0) &&
               (idIn1Abs == idNow[1]) ) skip = true;
        }
        if (skip) continue;
      }
      // Now decide which list of possible hadrons to use.
      // As we might have to use the special list for heavy quarks we
      // use the maximum of the absolute ids in case of two quarks and
      // check the maximum flavour in case of quark - diquark pair.
      int idUse;
      if ( (idIn1Abs < 10) && (idIn2Abs < 10) ) { // quark - quark
        idUse = ( (idIn1Abs > idIn2Abs) ? idIn1 : idIn2 );
      } else { // quark - diquark
        // Check if diquark contains a heavier flavour than the quark.
        bool useDiquark = false;
        for (int plus = 1; plus < 5; plus++)
          if (particleDataPtr->nQuarksInCode(idIn2Abs, idIn1Abs + plus) > 0)
            useDiquark = true;
        idUse = ( useDiquark ? idIn2 : idIn1 );
      }
      vector<double> possibleRatePrefacsNow = possibleRatePrefacs[idUse];
      vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idUse];
      // New list to fill.
      vector< pair<int,int> > possibleHadronsNew;
      vector<double> possibleRatePrefacsNew;
      // Now loop over possible hadrons and check if other (di)quark
      // in constituents matches idIn2.
      for (int iHad = 0; iHad < int(possibleHadronsNow.size()); iHad++) {
        // Get constituents.
        int hadronID = possibleHadronsNow[iHad].first;
        int iConst   = possibleHadronsNow[iHad].second;
        int ID1      = hadronConstIDs[hadronID][iConst].first;
        int ID2      = hadronConstIDs[hadronID][iConst].second;
        if ( ((ID1 == idIn1) && (ID2 == idIn2)) ||
             ((ID1 == idIn2) && (ID2 == idIn1)) ) {
          // Can take this combination.
          possibleHadronsNew.push_back(possibleHadronsNow[iHad]);
          possibleRatePrefacsNew.push_back(possibleRatePrefacsNow[iHad]);
        }
      }
      if (int(possibleHadronsNew.size()) < 1)
        loggerPtr->ERROR_MSG("no possible hadrons found for last two");
      // Save.
      possibleRatePrefacsLast[inPair] = possibleRatePrefacsNew;
      possibleHadronsLast[inPair]     = possibleHadronsNew;
    }
  }

  // Enhanced-rate prefactor for MPIs and/or nearby string pieces.
  closePacking     = flag("ClosePacking:doClosePacking");
  exponentMPI      = parm("ClosePacking:expMPI");
  exponentNSP      = parm("ClosePacking:expNSP");

  // Initialize winning parameters.
  hadronIDwin   = 0;
  quarkIDwin    = 0;
  hadronMassWin = -1.0;

}

//--------------------------------------------------------------------------

// Pick a hadron, based on generated pT value and initial (di)quark.
// Check all possible hadrons and calculate their relative suppression
// based on exp(-mThadron/T), possibly multiplied by spin counting, meson
// mixing or baryon weighting factors.
// First return value is hadron ID, second new (di)quark ID.

FlavContainer ThermalStringFlav::pick(FlavContainer& flavOld,
  double pT, double kappaModifier, bool) {

  // Determine close-packing scaling.
  double kappaRatio = 1. + enhancePT * kappaModifier;

  // Initial values for new flavour.
  FlavContainer flavNew;
  flavNew.rank = flavOld.rank + 1;

  int idIn        = flavOld.id;
  int idInAbs     = abs(idIn);
  double temprNow = temperature;
  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  if (idInAbs > 2) temprNow *= tempPreFactor;
  // Enhanced-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    temprNow     *= pow(max(1.0,kappaRatio), exponentNSP);
  }

  // Get the list of allowed hadrons and constituents for that
  // initial (di)quark. First parameter of pair is hadron ID, second
  // is nr of hadron constituents in the list.
  vector<double> possibleRatePrefacsNow      = possibleRatePrefacs[idIn];
  vector< pair<int,int> > possibleHadronsNow = possibleHadrons[idIn];
  int nPossHads = int(possibleHadronsNow.size());
  if (nPossHads < 1) {
    loggerPtr->ERROR_MSG("no possible hadrons found");
    return 0;
  }

  // Vector with hadron masses. Is -1.0 if m0 is use for calculating
  // the suppression rate and mSel if mSel is used.
  vector<double> possibleHadronMasses;

  // Calculate rates/suppression factors for given pT.
  vector<double> rates;
  double rateSum = 0.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    int hadronID = possibleHadronsNow[iHad].first;
    // Pick mass and calculate suppression factor.
    double mass  = particleDataPtr->mSel(hadronID);
    possibleHadronMasses.push_back(mass);
    double rate  = exp( -sqrt(pow2(pT)+pow2(mass))/temprNow );
    // Multiply rate with prefactor.
    rate *= possibleRatePrefacsNow[iHad];
    // Save rate and add to sum
    rates.push_back(rate);
    rateSum += rate;
  }
  // Normalize rates
  for (int iHad = 0; iHad < nPossHads; iHad++) rates[iHad] /= rateSum;

  // Get accumulated rates
  vector<double> accumRates;
  for (int iHad = 0; iHad < nPossHads; iHad++) accumRates.push_back(0);
  for (int iHad1 = 0; iHad1 < nPossHads; iHad1++)
    for (int iHad2 = 0; iHad2 <= iHad1; iHad2++)
      accumRates[iHad1] += rates[iHad2];

  // Random number to decide which hadron to pick
  double rand       = rndmPtr->flat();
  int hadronID      = 0;
  int iConst        = 0;
  double hadronMass = -1.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    if (rand <= accumRates[iHad]) {
      hadronID   = possibleHadronsNow[iHad].first;
      iConst     = possibleHadronsNow[iHad].second;
      hadronMass = possibleHadronMasses[iHad];
      break;
    }
  }

  // Get flavour of (di)quark to use next time.
  int idNext = 0;
  vector< pair<int,int> > constituentIDs = hadronConstIDs[hadronID];
  // Mesons
  if (particleDataPtr->isMeson(hadronID)) {
    int ID1 = constituentIDs[0].first;
    int ID2 = constituentIDs[0].second;
    // Special case for diagonal meson, flavour remains
    if (ID1 == -ID2) idNext = idIn;
    else idNext = (idIn == ID1 ? -ID2 : -ID1);
  }
  // Baryons
  else {
    int ID1 = constituentIDs[iConst].first;
    int ID2 = constituentIDs[iConst].second;
    if (ID1 == idIn) idNext = -ID2;
    if (ID2 == idIn) idNext = -ID1;
  }

  // Save new flavour and hadron.
  flavNew.id    = -idNext;  // id used to build hadron
  hadronIDwin   = hadronID;
  quarkIDwin    = idNext;   // id used in next step
  hadronMassWin = hadronMass;

  // Done.
  return flavNew;

}

//--------------------------------------------------------------------------

// Combine two flavours (including diquarks) to produce a hadron. Function
// called in case of combining the two remaining flavours into last hadron.

int ThermalStringFlav::combineLastThermal(FlavContainer& flav1,
  FlavContainer& flav2, double pT, double kappaModifier) {

  // Determine close-packing scaling.
  double kappaRatio = 1. + enhancePT * kappaModifier;

  // Decide randomly on whether to treat flav1 or flav2 as incoming.
  int idIn[2]    = { flav1.id, flav2.id };
  if (rndmPtr->flat() < 0.5) swap(idIn[0], idIn[1]);
  int idInNow[2] = { min(idIn[0],idIn[1]), max(idIn[0],idIn[1]) };

  int idInAbs     = abs(idIn[0]);
  double temprNow = temperature;
  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  if (idInAbs > 2) temprNow *= tempPreFactor;

  // Enhanced-rate prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow     *= pow(max(1.0,double(infoPtr->nMPI())), exponentMPI);
    temprNow     *= pow(max(1.0,kappaRatio), exponentNSP);
  }

  // Get the list of allowed hadrons and constituents for that combination
  // of (di)quarks. First parameter of pair is hadron ID, second
  // is nr of hadron constituents in the list.
  pair<int,int> inPr = pair<int,int>(idInNow[0], idInNow[1]);
  vector<double> possibleRatePrefacsNow      = possibleRatePrefacsLast[inPr];
  vector< pair<int,int> > possibleHadronsNow = possibleHadronsLast[inPr];
  int nPossHads = int(possibleHadronsNow.size());
  if (nPossHads < 1) {
    loggerPtr->ERROR_MSG("no possible hadrons found for last two");
    return 0;
  }

  // Vector with hadron masses. Is -1.0 if m0 is use for calculating
  // the suppression rate and mSel if mSel is used.
  vector<double> possibleHadronMasses;

  // Calculate rates/suppression factors for given pT.
  vector<double> rates;
  double rateSum = 0.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    int hadronID = possibleHadronsNow[iHad].first;
    // Pick mass and calculate suppression factor.
    double mass = particleDataPtr->mSel(hadronID);
    possibleHadronMasses.push_back(mass);
    double rate = exp( -sqrt(pow2(pT)+pow2(mass))/temprNow );
    // Multiply rate with prefactor.
    rate *= possibleRatePrefacsNow[iHad];
    // Save rate and add to sum
    rates.push_back(rate);
    rateSum += rate;
  }
  // Normalize rates
  for (int iHad = 0; iHad < nPossHads; iHad++) rates[iHad] /= rateSum;

  // Get accumulated rates
  vector<double> accumRates;
  for (int iHad = 0; iHad < nPossHads; iHad++) accumRates.push_back(0);
  for (int iHad1 = 0; iHad1 < nPossHads; iHad1++)
    for (int iHad2 = 0; iHad2 <= iHad1; iHad2++)
      accumRates[iHad1] += rates[iHad2];

  // Random number to decide which hadron to pick
  double rand       = rndmPtr->flat();
  int hadronID      = 0;
  double hadronMass = -1.0;
  for (int iHad = 0; iHad < nPossHads; iHad++) {
    if (rand <= accumRates[iHad]) {
      hadronID   = possibleHadronsNow[iHad].first;
      hadronMass = possibleHadronMasses[iHad];
      break;
    }
  }

  // Save hadron.
  hadronIDwin   = hadronID;
  hadronMassWin = hadronMass;

  // Done.
  return hadronIDwin;
}

//==========================================================================

// The ThermalStringPT class.

//--------------------------------------------------------------------------

// Initialize data members of the string pT selection.

void ThermalStringPT::init() {

  // Temperature for thermal model.
  temperature      = parm("StringPT:temperature");
  tempPreFactor    = parm("StringPT:tempPreFactor");

  // Upper estimate of thermal spectrum: fraction at x = pT_quark/T < 1.
  fracSmallX       = 0.6 / (0.6 + (1.2/0.9) * exp(-0.9));

  // Enhanced-width prefactor for MPIs and/or nearby string pieces.
  closePacking     = flag("ClosePacking:doClosePacking");
  enhancePT        = parm("ClosePacking:enhancePT");
  exponentMPI      = parm("ClosePacking:expMPI");
  exponentNSP      = parm("ClosePacking:expNSP");

}

//--------------------------------------------------------------------------

// Generate quark pT according to fitting functions, such that
// hadron pT is generated according to exp(-pT/T) d^2pT.

pair<double, double> ThermalStringPT::pxy(int idIn, double kappaModifier) {

  // Temperature increase to work against asymmetry. Apply for
  // s/c/b and diquarks.
  double temprNow = temperature;
  if (abs(idIn) > 2) temprNow *= tempPreFactor;

  // Enhanced-width prefactor for MPIs and/or nearby string pieces.
  if (closePacking) {
    temprNow *= pow(max(1.0, double(infoPtr->nMPI())), exponentMPI);
    temprNow *= pow(max(1.0, kappaModifier), exponentNSP);
  }

  // Pick x = pT_quark/T according to K_{1/4}(x)/x^{1/4} * x dx.
  double xrand, approx, wanted;
  do {
    xrand = (rndmPtr->flat() < fracSmallX) ? rndmPtr->flat()
          : 1. - log(rndmPtr->flat()) / 0.9;
    approx = (xrand < 1.) ? 0.6 : 1.2 * exp(-0.9 * xrand);
    wanted = BesselK14(xrand) * pow( xrand, 0.75);
  } while (rndmPtr->flat() * approx > wanted);

  // Find pT_quark. Random number to decide on angle.
  double pTquark = xrand * temprNow;
  double phi     = 2.0 * M_PI * rndmPtr->flat();

  // Done.
  return pair<double, double>( pTquark * cos(phi), pTquark * sin(phi) );

}

//--------------------------------------------------------------------------

// Evaluate Bessel function K_{1/4}(x).
// Use power series for x < 2.5 and asymptotic expansion for x > 2.5.
// Number of terms picked to have accuracy better than 1 per mille.
// Based on M. Abramowitz and I.A. Stegun, eqs. 9.6.2, 9.6.10, 9.7.2.

double ThermalStringPT::BesselK14(double x) {

  // Power series expansion of K_{1/4} : k = 0 term.
  if (x < 2.5) {
    double xRat  = 0.25 * x * x;
    double prodP = pow( 0.5 * x, -0.25) / 1.2254167024;
    double prodN = pow( 0.5 * x,  0.25) / 0.9064024771;
    double sum   = prodP - prodN;

    // Power series expansion of K_{1/4} : m > 0 terms.
    for (int k = 1; k < 6; ++k) {
      prodP *= xRat / (k * (k - 0.25));
      prodN *= xRat / (k * (k + 0.25));
      sum   += prodP - prodN;
    }
    sum *= M_PI * sqrt(0.5);
    return sum;

  // Asymptotic expansion of K_{1/4}.
  } else {
    double asym  = sqrt(M_PI * 0.5 / x) * exp(-x);
    double term1 = -         0.75 / ( 8. * x);
    double term2 = -term1 *  8.75 / (16. * x);
    double term3 = -term2 * 24.75 / (24. * x);
    double term4 = -term3 * 48.75 / (32. * x);
    asym *= 1. + term1 + term2 + term3 + term4;
    return asym;
  }
}

//==========================================================================

// The ThermalFragmentation class.

//--------------------------------------------------------------------------

// Constructor.

ThermalFragmentation::ThermalFragmentation() {

  thermalFlavSelPtr = new ThermalStringFlav();
  thermalPTSelPtr   = new ThermalStringPT();
  zSelPtr           = new StringZ();
  stringFragPtr     = new StringFragmentation();
  ministringFragPtr = new MiniStringFragmentation();

}

//--------------------------------------------------------------------------

// Destructor.

ThermalFragmentation::~ThermalFragmentation() {

  delete stringFragPtr;
  delete ministringFragPtr;
  delete thermalFlavSelPtr;
  delete thermalPTSelPtr;
  delete zSelPtr;

}

//--------------------------------------------------------------------------

// Initialize and save pointers.

bool ThermalFragmentation::init(StringFlav*, StringPT*, StringZ*,
  FragModPtr fragModPtrIn) {

  // Register the string and ministring fragmentation models.
  registerSubObject(*thermalFlavSelPtr);
  registerSubObject(*thermalPTSelPtr);
  registerSubObject(*zSelPtr);
  registerSubObject(*stringFragPtr);
  registerSubObject(*ministringFragPtr);

  // Initialize the ThermalStringFlav, ThermalStringPT and zSelPtr objects.
  thermalFlavSelPtr->init();
  thermalPTSelPtr->init();
  zSelPtr->init();

  // Set the pointers.
  stringFragPtr->setFlavBeforePT( false);
  stringFragPtr->init(thermalFlavSelPtr, thermalPTSelPtr, zSelPtr,
    fragModPtrIn);
  ministringFragPtr->init(thermalFlavSelPtr, thermalPTSelPtr, zSelPtr,
    fragModPtrIn);

  // Boundary mass between string and ministring handling.
  mStringMin = parm("HadronLevel:mStringMin");

  // Try ministring fragmentation also if normal fails.
  tryMiniAfterFailedFrag = flag("MiniStringFragmentation:tryAfterFailedFrag");

  // Return successful.
  return true;

}

//--------------------------------------------------------------------------

// Do the fragmentation: driver routine.

bool ThermalFragmentation::fragment(int iSub, ColConfig& colConfig,
  Event& event, bool isDiff, bool) {

  // String fragmentation of each colour singlet (sub)system.
  // If fails optionally try ministring fragmentation.
  if (iSub == -1) return true;
  if (colConfig[iSub].massExcess > mStringMin) {
    if (!stringFragPtr->fragment( iSub, colConfig, event)) {
      if (!tryMiniAfterFailedFrag) return false;
      loggerPtr->ERROR_MSG("string fragmentation failed, "
        "trying ministring fragmetation instead");
      if (!ministringFragPtr->fragment(iSub, colConfig, event, isDiff)) {
        loggerPtr->ERROR_MSG("also ministring fragmentation failed "
          "after failed normal fragmentation");
        return false;
      }
    }

  // Low-mass string treated separately.
  } else {
    if (!ministringFragPtr->fragment( iSub, colConfig, event, isDiff)) {
      loggerPtr->ERROR_MSG("ministring fragmentation failed");
      return false;
    }
  }

  // Return successful.
  return true;

}

//==========================================================================

} // end namespace Pythia8
