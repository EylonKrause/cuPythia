// main424.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Marius Utheim
//          Leif Lonnblad <leif.lonnblad@fysik.lu.se>

// Keywords: hadron-ion collisions, cosmic rays, optimization

// The purpose of this example is to generate an initialization file that
// can be used to speed up initialization in hadron-hadron or hadron-ion runs.
// By default, it produces data for energies from 10 to 10^7 GeV.
// All hadron-nucleon and hadron-ion interactions are possible.
// In other main programs, after initializing with the main424.cmnd file,
// it is possible to change energy and beam types on a per-event basis.
// Warning: the generation takes a long time, of the order of 3 hours,
// whereof most of it is spent on the nuclear geometry initialization.
// But, of course, once done you will save time in subsequent runs.

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(int argc, char ** argv) {

  int option = argc > 1? atoi(argv[1]): 0;
  // option = 0: generate combined MPI and Angantyr file.
  // option = 1: generate only MPI file.
  // option = 2: generate only Angantyr file (assumes MPI file exists).
  string initMPIcmnd = "include = setups/InitDefaultMPI.cmnd";
  string filename = "main424";
  if ( option > 0 ) filename += string("-") + argv[1];
  filename += ".cmnd";

  // The maximum energy.
  double eCMMax = 1e7;

  // The number of points in the (logarithmic) interpoaltion grid for
  // Angantyr (can be lowered to get faster results):
  int nGridHI   = 25;

  // The number of integrations points for the MPI interpolation.
  int nIntPts   = 1000000;

  // Use Angantyr even when initializing with pp.
  Pythia pythia;

  // Variable energy parameters.
  pythia.readString("Beams:allowVariableEnergy = on");
  pythia.settings.parm("Beams:eCM", eCMMax);
  pythia.readString("SoftQCD:all = on");

  // Variable beam type parameters.
  pythia.readString("Beams:allowIDASwitch = on");

  if ( option != 1 ) {
    // Setup Angantyr for initialisation.
    pythia.readString("HeavyIon:mode = 2");
    pythia.settings.parm("HeavyIon:varECMMax", eCMMax);
    pythia.settings.parm("HeavyIon:varECMMin", 10.);
    pythia.settings.mode("HeavyIon:varECMSigFitNPts", nGridHI);

    // For subcollision models other than the default, use 1 here:
    pythia.settings.mode("HeavyIon:varECMStepwiseEvolve", 2);

    // Set reasonable starting parameters for the maximum energy. And
    // reasonable target errors for semi-inclusive xsecs (noting that
    // double diffrection is almost impossible to fit).
    pythia.readString("HeavyIon:SigFitDefPar = 4.73,32.65,0.69");
    pythia.readString("HeavyIon:SigFitErr ="
       "0.02,0.02,0.0,0.02,0.02,0.0,0.02,0.2");
    pythia.readString("HeavyIon:SigFitNGen = 20");
    pythia.readString("HeavyIon:SigFitMaxChi2 = 1.5");
    pythia.readString("HeavyIon:SigFitMaxChi2Max = 4");

    // If you set reuseInit > 0, you also need to specify file names
    // here. Otherwise the files are produced explicitly below.
    pythia.readString("HeavyIon:SasdMpiReuseInit = -1");
    pythia.readString("HeavyIon:SasdMpiInitFile = ");
    pythia.readString("HeavyIon:SigFitReuseInit = -1");
    pythia.readString("HeavyIon:SigFitInitFile = ");
  }

  if ( option == 2 ) {
    // Read in previous MPI file.
    pythia.readString(initMPIcmnd);
  } else {
    // Increased statistics for MPI initialization gives smaller errors.
    pythia.settings.mode("MultipartonInteractions:nSample", nIntPts);

    // If you set reuseInit > 0, you also need to specify file names
    // here. Otherwise the files are produced explicitly below.
    pythia.readString("MultipartonInteractions:reuseInit = -1");
    pythia.readString("MultipartonInteractions:initFile = ");
  }

  // Also save and print the maximum energy for easy checks.
  pythia.settings.parm("Beams:eCMMaxMPI", eCMMax);

  // Specify where to save. If you set reuseInit = 2, the old files
  // will be replaced/overwritten if they already exist. -1 means that
  // the initialization is stored in the settings, and can be written
  // out.

  if (!pythia.init()) {
    cout << " Pythia failed to initialize." << endl
         << " Please validate your settings in main424.cmnd." << endl;
    return -1;
  }

  // Note that writeFile will write out all changed settings, not only
  // the ones containing the initialisations, so some manual clean-up
  // is needed to avoid this.
  pythia.readString("Beams:eCM = DEFAULT");
  pythia.readString("HeavyIon:mode = DEFAULT");
  pythia.readString("Next:numberCount = DEFAULT");
  pythia.readString("Next:numberShowEvent = DEFAULT");
  pythia.readString("Next:numberShowInfo = DEFAULT");
  pythia.readString("Next:numberShowLHA = DEFAULT");
  pythia.readString("Next:numberShowProcess = DEFAULT");
  pythia.readString("PartonLevel:all = DEFAULT");
  pythia.readString("ProcessLevel:all = DEFAULT");
  pythia.readString("SoftQCD:all = DEFAULT");
  pythia.readString("Tune:ee = DEFAULT");
  pythia.readString("Tune:pp = DEFAULT");
  pythia.readString("StringZ:useOldAExtra = DEFAULT");
  pythia.readString("MultipartonInteractions:initFile = DEFAULT");
  pythia.readString("HeavyIon:SigFitInitFile = DEFAULT");
  pythia.readString("HeavyIon:SasdMpiInitFile = DEFAULT");
  pythia.readString("HeavyIon:SigFitDefPar = DEFAULT");
  pythia.readString("HeavyIon:SigFitErr = DEFAULT");
  pythia.readString("HeavyIon:SigFitNGen = DEFAULT");
  pythia.readString("HeavyIon:SigFitMaxChi2 = DEFAULT");
  pythia.readString("HeavyIon:SigFitMaxChi2Max = DEFAULT");
  pythia.readString("HeavyIon:varECMMax = DEFAULT");
  pythia.readString("HeavyIon:mode = DEFAULT");

  if ( option == 2 ) {
    pythia.readString("Init:reuseMPIiDiffSys0 = DEFAULT");
    pythia.readString("Init:reuseMPIiDiffSys1 = DEFAULT");
    pythia.readString("Init:reuseMPIiDiffSys2 = DEFAULT");
    pythia.readString("Init:reuseMPIiDiffSys3 = DEFAULT");
  }

  ofstream os(filename.c_str());
  if ( option == 2 ) os << initMPIcmnd << endl;
  pythia.settings.writeFile(os);

  // After initializing, we're done.

  return 0;
}
