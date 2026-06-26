// HeavyIons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HeavyIons.h header) for the
// heavy ion classes classes, and some related global functions.

#include "Pythia8/BeamShape.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8/HINucleusModel.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// The abstract HeavyIons class

//--------------------------------------------------------------------------

// Before doing anything, Pythia should add special heavy ion versions
// for some groups of settings, to be used in the generation of
// secondary absorptive scatterings in case Angantyr:SASDmode > 0.

void HeavyIons::addSpecialSettings(Settings & settings) {
  setupSpecials(settings, "Diffraction:");
  setupSpecials(settings, "MultipartonInteractions:");
  setupSpecials(settings, "PDF:");
  setupSpecials(settings, "SigmaDiffractive:");
  setupSpecials(settings, "BeamRemnants:");
  settings.addMode("Angantyr:testMode", 0, false, false, 0, 0);
}

//--------------------------------------------------------------------------

// Duplicate all settings of the form "<match>..." to "HI<...>" variants,
// for use by HeavyIon sub-Pythia objects. Needs updating if new Settings
// categories are introduced.

void HeavyIons::setupSpecials(Settings & settings, string match) {
  map<string,Flag> flags = settings.getFlagMap(match);
  for ( map<string,Flag>::iterator it = flags.begin();
        it != flags.end(); ++it )
    settings.addFlag("HI" + it->second.name, it->second.valDefault);
  map<string,Mode> modes = settings.getModeMap(match);
  for ( map<string,Mode>::iterator it = modes.begin();
        it != modes.end(); ++it )
    settings.addMode("HI" + it->second.name, it->second.valDefault,
                     it->second.hasMin, it->second.hasMax,
                     it->second.valMin, it->second.valMax, it->second.optOnly);
  map<string,Parm> parms = settings.getParmMap(match);
  for ( map<string,Parm>::iterator it = parms.begin();
        it != parms.end(); ++it )
    settings.addParm("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,Word> words = settings.getWordMap(match);
  for ( map<string,Word>::iterator it = words.begin();
        it != words.end(); ++it )
    settings.addWord("HI" + it->second.name, it->second.valDefault);
  map<string,FVec> fvecs = settings.getFVecMap(match);
  for ( map<string, FVec>::iterator it = fvecs.begin();
        it != fvecs.end(); ++it )
    settings.addFVec("HI" + it->second.name, it->second.valDefault);
  map<string,MVec> mvecs = settings.getMVecMap(match);
  for ( map<string,MVec>::iterator it = mvecs.begin();
        it != mvecs.end(); ++it )
    settings.addMVec("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,PVec> pvecs = settings.getPVecMap(match);
  for ( map<string,PVec>::iterator it = pvecs.begin();
        it != pvecs.end(); ++it )
    settings.addPVec("HI" + it->second.name, it->second.valDefault,
                 it->second.hasMin, it->second.hasMax,
                 it->second.valMin, it->second.valMax);
  map<string,WVec> wvecs = settings.getWVecMap(match);
  for ( map<string,WVec>::iterator it = wvecs.begin();
        it != wvecs.end(); ++it )
    settings.addWVec("HI" + it->second.name, it->second.valDefault);
}

void HeavyIons::setupSpecials(Pythia & p, string match) {
  Settings & opts = p.settings;
  map<string, Flag> flags = opts.getFlagMap(match);
  for ( map<string, Flag>::iterator it = flags.begin();
        it != flags.end(); ++it )
    opts.flag(it->second.name.substr(2), it->second.valNow, true);
  map<string, Mode> modes = opts.getModeMap(match);
  for ( map<string, Mode>::iterator it = modes.begin();
        it != modes.end(); ++it )
    opts.mode(it->second.name.substr(2), it->second.valNow, true);
  map<string, Parm> parms = opts.getParmMap(match);
  for ( map<string, Parm>::iterator it = parms.begin();
        it != parms.end(); ++it )
    opts.parm(it->second.name.substr(2), it->second.valNow, true);
  map<string, Word> words = opts.getWordMap(match);
  for ( map<string, Word>::iterator it = words.begin();
       it != words.end(); ++it )
    opts.word(it->second.name.substr(2), it->second.valNow, true);
  map<string, FVec> fvecs = opts.getFVecMap(match);
  for ( map<string, FVec>::iterator it = fvecs.begin();
        it != fvecs.end(); ++it )
    opts.fvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, MVec> mvecs = opts.getMVecMap(match);
  for ( map<string, MVec>::iterator it = mvecs.begin();
        it != mvecs.end(); ++it )
    opts.mvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, PVec> pvecs = opts.getPVecMap(match);
  for ( map<string, PVec>::iterator it = pvecs.begin();
        it != pvecs.end(); ++it )
    opts.pvec(it->second.name.substr(2), it->second.valNow, true);
  map<string, WVec> wvecs = opts.getWVecMap(match);
  for ( map<string, WVec>::iterator it = wvecs.begin();
        it != wvecs.end(); ++it )
    opts.wvec(it->second.name.substr(2), it->second.valNow, true);
}

//--------------------------------------------------------------------------

// Clear SoftQCD flags for Pythia subobjects, but return processes
// explicitly turned off.

set<int> HeavyIons::clearSoftQCDFlags(Settings& settings) {

  // Map per-channel SoftQCD and LowEnergyQCD flags to internal process codes.
  const static map<string,int> softMap = {
    {"nonDiffractive",      101},
    {"elastic",             102},
    {"singleDiffractiveXB", 103},
    {"singleDiffractiveAX", 104},
    {"doubleDiffractive",   105},
    {"centralDiffractive",  106}
  };

  const static map<string,int> lowMap = {
    {"nonDiffractive",      151},
    {"elastic",             152},
    {"singleDiffractiveXB", 153},
    {"singleDiffractiveAX", 154},
    {"doubleDiffractive",   155},
    {"excitation",          157},
    {"annihilation",        158},
    {"resonant",            159}
  };

  // Build veto list for SoftQCD processes.
  set<int> vetoSoft;

  if (settings.flag("SoftQCD:all")) {
    // User explicitly allows all SoftQCD. Nothing is vetoed.
  } else if (settings.flag("SoftQCD:inelastic")) {
    // Special case: "inelastic" means elastic-only is off.
    vetoSoft.insert(102);
  } else {
    // Start from "everything off" and re-enable user-requested channels.
    for (const auto& kv : softMap)
      vetoSoft.insert(kv.second);

    int count = 0;
    // Remove channels explicitly turned on.
    for (const auto& kv : softMap)
      if (settings.flag("SoftQCD:" + kv.first)) {
        vetoSoft.erase(kv.second);
        ++count;
      }
    // For historic reasons we allow the user to not specify any
    // SoftQCD processes, in which case they get all of them.
    if ( !count ) vetoSoft.clear();
  }
  // Repeat for LowEnergyQCD processes.
  set<int> vetoLow;

  if (settings.flag("LowEnergyQCD:all")) {
    // All allowed; nothing vetoed.
  } else {
    // Start from full set.
    for (const auto& kv : lowMap)
      vetoLow.insert(kv.second);

    // Remove explicitly enabled channels.
    for (const auto& kv : lowMap)
      if (settings.flag("LowEnergyQCD:" + kv.first))
        vetoLow.erase(kv.second);
  }
  // Reset all per-channel flags.
  // SoftQCD includes one special aggregate flag (inelastic).
  settings.flag("SoftQCD:inelastic", false);
  for (const auto& kv : softMap)
    settings.flag("SoftQCD:" + kv.first, false);

  // LowEnergyQCD: reset all per-channel flags from lowMap keys.
  for (const auto& kv : lowMap)
    settings.flag("LowEnergyQCD:" + kv.first, false);

  // Merge vetoes and return.
  vetoSoft.insert(vetoLow.begin(), vetoLow.end());
  return vetoSoft;

}

//--------------------------------------------------------------------------

// Reset all process level settings in the given Pythia object. NOTE
// must be expanded if new process groups are included in Pythia.

void HeavyIons::clearProcessLevel(Pythia & pyt) {
  string path = pyt.settings.word("xmlPath");
  pyt.settings.mode("Tune:ee", 0);
  pyt.settings.mode("Tune:pp", 0);
  pyt.settings.init(path + "QCDSoftProcesses.xml", true);
  pyt.settings.init(path + "QCDHardProcesses.xml", true);
  pyt.settings.init(path + "ElectroweakProcesses.xml", true);
  pyt.settings.init(path + "OniaProcesses.xml", true);
  pyt.settings.init(path + "TopProcesses.xml", true);
  pyt.settings.init(path + "FourthGenerationProcesses.xml", true);
  pyt.settings.init(path + "HiggsProcesses.xml", true);
  pyt.settings.init(path + "SUSYProcesses.xml", true);
  pyt.settings.init(path + "NewGaugeBosonProcesses.xml", true);
  pyt.settings.init(path + "LeftRightSymmetryProcesses.xml", true);
  pyt.settings.init(path + "LeptoquarkProcesses.xml", true);
  pyt.settings.init(path + "CompositenessProcesses.xml", true);
  pyt.settings.init(path + "HiddenValleyProcesses.xml", true);
  pyt.settings.init(path + "ExtraDimensionalProcesses.xml", true);
  pyt.settings.init(path + "DarkMatterProcesses.xml", true);
  pyt.settings.init(path + "SecondHardProcess.xml", true);
  pyt.settings.init(path + "PhaseSpaceCuts.xml", true);
  // NOTE! if new processes are added in separate xml files these have
  // to be added here.
}

//--------------------------------------------------------------------------

// Update the Info object in the main Pythia object.

void HeavyIons::updateInfo() {
  *infoPtr = hiInfo.primInfo;
  infoPtr->particleDataPtr = particleDataPtr;
  infoPtr->hiInfo = &hiInfo;
  infoPtr->weightContainerPtr->setWeightNominal(hiInfo.weight());
  infoPtr->sigmaReset();
  double norm = 1.0/double(hiInfo.NSave);
  int Nall = 0;
  double wall = 0.0;
  double w2all = 0.0;
  for ( map<int,int>::iterator ip = hiInfo.NPrim.begin();
        ip != hiInfo.NPrim.end(); ++ip ) {
    int N = ip->second;
    if ( !N ) continue;
    int pc = ip->first;
    double w = hiInfo.sumPrimW[pc]*FMSQ2MB;
    double w2 = hiInfo.sumPrimW2[pc]*pow2(FMSQ2MB);
    infoPtr->setSigma(pc, hiInfo.NamePrim[pc], N, N, N,
                      w*norm, sqrt(w2*norm)/N, w * MB2FMSQ);
    Nall += N;
    wall += w;
    w2all += w2;
  }
  infoPtr->setSigma(0, "sum", hiInfo.NSave, Nall, Nall,
                    wall*norm, sqrt(w2all*norm)/Nall, wall * MB2FMSQ);
}

//--------------------------------------------------------------------------

// Print out statistics from a HeavyIons run.

void HeavyIons::stat() {
  bool showPrL = flag("Stat:showProcessLevel");
  bool showErr = flag("Stat:showErrors");
  bool reset   = flag("Stat:reset");
  Info & in = *infoPtr;
  // Header.
  if ( showPrL ) {
    cout << "\n *-----  HeavyIon Event and Cross Section Statistics  ------"
         << "-------------------------------------------------------*\n"
         << " |                                                            "
         << "                                                     |\n"
         << " | Primary NN sub-collision subprocess           Code |       "
         << "     Number of events       |      sigma +- delta    |\n"
         << " |                                                    |       "
         << "Tried   Selected   Accepted |     (estimated) (mb)   |\n"
         << " |                                                    |       "
         << "                            |                        |\n"
         << " |------------------------------------------------------------"
         << "-----------------------------------------------------|\n"
         << " |                                                    |       "
         << "                            |                        |\n";

    vector<int> pc = in.codesHard();
    bool caveat = false;
    for ( int i = 0, N = pc.size(); i < N; ++i ) {
      string pname = in.nameProc(pc[i]);
      if ( pc[i] == 102 ) {
        pname += " (*)";
        caveat = true;
      }
      cout << " | " << left << setw(45) << pname
           << right << setw(5) << pc[i] << " | "
           << setw(11) << in.nTried(pc[i]) << " "
           << setw(10) << in.nSelected(pc[i]) << " "
           << setw(10) << in.nAccepted(pc[i]) << " | "
           << scientific << setprecision(3)
           << setw(11) << in.sigmaGen(pc[i])
           << setw(11) << in.sigmaErr(pc[i]) << " |\n";
    }
    if ( pc.empty() ) in.setSigma(0, "sum", hiInfo.NSave, 0, 0, 0.0, 0.0, 0.0);

    cout << " | " << left << setw(50) << "sum" << right << " | " << setw(11)
         << in.nTried(0) << " " << setw(10) << in.nSelected(0) << " "
         << setw(10) << in.nAccepted(0) << " | " << scientific
         << setprecision(3) << setw(11)
         << in.sigmaGen(0) << setw(11) << in.sigmaErr(0) << " |\n"
         << " |                                                    |       "
         << "                            |                        |\n";
    if ( caveat )
      cout << " | (*) Note: elastic events are not correctly treated |       "
           << "                            |                        |\n";
    if ( hasGlauberCalculation() ) {
      cout << " |------------------------------------------------------------"
           << "-----------------------------------------------------|\n"
           << " |                                                            "
           << "                            |                        |\n";
      string line = "Semi-inclusive " + particleDataPtr->name(idProj) + " on "
        + particleDataPtr->name(idTarg) +
        " cross sections from the Glauber calculation:";
      cout << " | " << left << setw(86) << line
           << " |                        |\n"
           << " |                                                            "
           << "                            |                        |\n";
      cout << " | " << left << setw(86)
           << "Total" << " | "
           << right << scientific << setprecision(3)
           << setw(11) << hiInfo.glauberTot()
           << setw(11) << hiInfo.glauberTotErr() << " |\n";
    }
    if ( hasGlauberCalculation() > 1 ) {
      cout << " | " << left << setw(86)
           << "Non-Diffractive" << " | "
           << right << scientific << setprecision(3)
           << setw(11) << hiInfo.glauberND()
           << setw(11) << hiInfo.glauberNDErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Total inelastic" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberINEL() << setw(11)
           << hiInfo.glauberINELErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Elastic" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberEL() << setw(11)
           << hiInfo.glauberELErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Diffractive target excitation" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberDiffT() << setw(11)
           << hiInfo.glauberDiffTErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Diffractive projectile excitation" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberDiffP() << setw(11)
           << hiInfo.glauberDiffPErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Double diffractive excitation" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberDDiff() << setw(11)
           << hiInfo.glauberDDiffErr() << " |\n";
      cout << " | " << left << setw(86)
           << "Elastic b-slope (GeV^-2)" << " | "
           << right << scientific << setprecision(3) << setw(11)
           << hiInfo.glauberBSlope() << setw(11)
           << hiInfo.glauberBSlopeErr() << " |\n";
    }
    // Listing finished.
    cout << " |                                                            "
         << "                            |                        |\n"
         << " *-----  End HeavyIon Event and Cross Section Statistics -----"
         << "-----------------------------------------------------*" << endl;
  }
  if ( reset ) hiInfo = HIInfo();
  if ( showErr ) {
    for ( int i = 1, np = pythia.size(); i < np; ++i )
      loggerPtr->errorCombine(pythia[i]->logger, "(" + pythiaNames[i] + ")");
    loggerPtr->errorStatistics();
  }
  if ( reset ) loggerPtr->errorReset();

}

//--------------------------------------------------------------------------

// Check the settings and return false of there are no heavy ion beams.

bool HeavyIons::isHeavyIon(Settings & settings) {
  int idProj = settings.mode("Beams:idA");
  int idTarg = settings.mode("Beams:idB");
  return ( abs(idProj/100000000) == 10 ||abs(idTarg/100000000) == 10 );
}

//==========================================================================

// Angantyr is the main HeavyIons model in Pythia.

//--------------------------------------------------------------------------

// Constructor.

Angantyr::Angantyr(Pythia & mainPythiaIn)
  : HeavyIons(mainPythiaIn) {
  selectMB = make_shared<ProcessSelectorHook>();
  selectSASD = make_shared<ProcessSelectorHook>();
  pythia.resize(ALL);
  info.resize(ALL);
  pythiaNames.resize(ALL);
  pythiaNames[HADRON] = "HADRON";
  pythiaNames[MBIAS] = "MBIAS";
  pythiaNames[SASD] = "SASD";
  pythiaNames[SIGPP] = "SIGPP";
  pythiaNames[SIGPN] = "SIGPN";
  pythiaNames[SIGNP] = "SIGNP";
  pythiaNames[SIGNN] = "SIGNN";

}

//--------------------------------------------------------------------------

// Destructor deleting model objects that are not provided from the
// outside (via HIUserHooks).

Angantyr::~Angantyr() {
  for ( int i = MBIAS; i < ALL; ++i ) if ( pythia[i] ) delete pythia[i];
}

//--------------------------------------------------------------------------

// Add a HIUserHooks object to customise the Angantyr model.

bool Angantyr::setUserHooksPtr(PythiaObject sel, shared_ptr<UserHooks> uhook) {
  for ( int i = HADRON; i < ALL; ++i )
    if ( ( i == sel || ALL == sel ) && !pythia[i]->setUserHooksPtr(uhook) )
      return false;
  return true;
}

//--------------------------------------------------------------------------

// Figure out what beams the user wants.

void Angantyr::setBeamKinematics(int idA, int idB) {
  // We will use the MBIAS BeamSetup object to figure out what is
  // happening. Whatever we do here will be overridden when we do the
  // proper init().
  pythia[MBIAS]->settings.mode("Beams:idA", idA);
  pythia[MBIAS]->settings.mode("Beams:idB", idB);
  beamSetupPtr->mA = particleDataPtr->m0(idA);
  beamSetupPtr->mB = particleDataPtr->m0(idB);
  if ( idProj != idA ) {
    int A = (idProj/10)%1000;
    beamSetupPtr->mA = particleDataPtr->m0(idProj)/A;
  }
  if ( idTarg != idB ) {
    int A = (idTarg/10)%1000;
    beamSetupPtr->mB = particleDataPtr->m0(idTarg)/A;
  }
  beamSetupPtr->initFrame();
  unifyFrames();
}

//--------------------------------------------------------------------------

// Update the medium cross section overestimates for Cascade mode.

void Angantyr::updateMedium() {
  int idBSave = targPtr->id();
  hiInfo.mediumXSecsSave.clear();
  for ( int idMedium : cascadeMediumIons ) {
    targPtr->setParticle(idMedium);
    bGenPtr->updateWidth();
    hiInfo.mediumXSecsSave.push_back(bGenPtr->xSecScale()*FMSQ2MB);
  }
  targPtr->setParticle(idBSave);
}

//--------------------------------------------------------------------------

// Switch to new beam particle identities.

bool Angantyr::setBeamIDs(int idAIn, int idBIn) {

  if ( idAIn == projPtr->id() && ( idBIn == 0 || idBIn == targPtr->id() ) )
    return true;

  if ( idBIn == 0 ) idBIn = targPtr->id();

  // Reset the statistics.
  hiInfo.glauberReset();

  // Set the projectile and target IDs.
  projPtr->setParticle(idAIn);
  targPtr->setParticle(idBIn);

  // Set the beam IDs in minimum bias.
  if (!pythia[MBIAS]->setBeamIDs(projPtr->idN(), targPtr->idN()))
    return false;
  if (sabsMode >= 0 &&
      !pythia[SASD]->setBeamIDs(projPtr->idN(), targPtr->idN()))
    return false;

  // Set masses and IDs.
  beamSetupPtr->mA = projPtr->mN();
  beamSetupPtr->mB = targPtr->mN();
  beamSetupPtr->idA = idAIn;
  beamSetupPtr->idB = idBIn;
  unifyFrames();

  // Beam masses may have changed so we reset kinematics as well.
  setKinematicsCM();

  // Calculate the total cross-section.
  if ( !doLowEnergyNow )
    sigTotNN.calc(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM);
  else if ( !lowEnergyCollPtr->hasXSec() )
    return false;

  // Update Subcollision and impact parameter handlers.
  if ( !collPtr->setIDA(beamSetupPtr->represent(projPtr->idN())) )
    return false;

  // If cascade mode update cross sections for all nuclei in the
  // medium.
  if ( cascadeMode ) updateMedium();

  bGenPtr->updateWidth();

  idProj = idAIn;
  idTarg = idBIn;

  return true;
}

//--------------------------------------------------------------------------

// Create an EventInfo object connected to a SubCollision from the
// last event generated by the given PythiaObject.

EventInfo Angantyr::mkEventInfo(Pythia & pyt, Info & infoIn,
                                const SubCollision * coll) {
  EventInfo ei;
  ei.coll = coll;
  ei.event = pyt.event;
  ei.info = infoIn;
  ei.code =  pyt.info.code();
  ei.ordering = ( ( HIHooksPtr && HIHooksPtr->hasEventOrdering() )?
                  HIHooksPtr->eventOrdering(ei.event, infoIn):
                  pyt.info.bMPI() );
  if ( coll ) {
    ei.projs[coll->proj] = make_pair(1, ei.event.size());
    ei.targs[coll->targ] = make_pair(2, ei.event.size());
  }

  ei.ok = streamline(ei);

  return ei;
}

//--------------------------------------------------------------------------

// Display the Angantyr banner.

void Angantyr::banner() const {

  string colOut = "              ";
  string cols = particleDataPtr->name(idProj)+" on "+
    particleDataPtr->name(idTarg);
  colOut.replace(colOut.begin(), colOut.begin() + cols.size(), cols);

  cout << " *----------------------  Initializing Angantyr  ----------------"
       << "------*\n"
       << " |                    We collide: " + colOut + "                 "
       << "      |\n"
       << " |                                                               "
       << "      |\n"
       << " |                    Below follows initialization               "
       << "      |\n"
       << " |                    of sub-collisions modelling.               "
       << "      |\n"
       << " |                                                               "
       << "      |\n"
       << " |                   //>________________________________         "
       << "      |\n"
       << " |          [########[]_________________________________>        "
       << "      |\n"
       << " |                   \\\\>                                       "
       << "        |\n"
       << " *-------------------------------------------------------------"
       << "--------*" << endl;

}

//--------------------------------------------------------------------------

// Initialise Angantyr. Called from within Pythia::init().

bool Angantyr::init() {

  // Read settings.
  idProj         = mode("Beams:idA");
  idTarg         = mode("Beams:idB");
  doSDTest       = flag("Angantyr:SDTest");
  glauberOnly    = flag("Angantyr:GlauberOnly");
  recoilerMode   = mode("Angantyr:SDRecoil");
  bMode          = mode("Angantyr:impactMode");
  doVarECM       = flag("Beams:allowVariableEnergy");
  eCMlow         = parm("HeavyIon:eCMLowEnergy");
  doHadronLevel  = flag("HadronLevel:all");
  sabsMode       = mode("Angantyr:SASDmode");
  sabsEps        = parm("Angantyr:epsilonSABS");
  sabsMinMX      = parm("Angantyr:minMxSABS");
  sabsCutMX      = parm("Angantyr:cutMxSABS");
  allowIDAswitch = flag("Beams:allowIDAswitch");

  if ( ( cascadeMode = flag("Angantyr:cascadeMode") ) ) {
    if ( mode("Angantyr:CollisionModel") != 6 ||
         !flag("HeavyIon:forceUnitWeight")) {
      loggerPtr->ABORT_MSG("Angantyr:cascadeMode=on requires "
                           "Angantyr:CollisionModel=6 and "
                           "HeavyIon:forceUnitWeight=on");
      return false;
    }
    cascadeMediumIons = mvec("Angantyr:cascadeMediumIons");
  }

  int idProjP = idProj;
  int idProjN = 0;
  int idTargP = idTarg;
  int idTargN = 0;
  bool isHIProj = ( abs(idProj/100000000) == 10 );
  bool isHITarg = ( abs(idTarg/100000000) == 10 );
  bool isHI = isHIProj || isHITarg || mode("HeavyIon:mode") > 1;
  if ( isHIProj ) {
    idProjP = idProj > 0? 2212: -2212;
    idProjN = idProj > 0? 2112: -2112;
  }
  if ( isHITarg ) {
    idTargP = idTarg > 0? 2212: -2212;
    idTargN = idTarg > 0? 2112: -2112;
  }
  if ( mode("HeavyIon:mode") == 1 && !isHI ) {
    loggerPtr->ABORT_MSG("no heavy ions requested");
    return false;
  }

  bool print = flag("HeavyIon:showInit") && !settingsPtr->flag("Print:quiet");
  if ( print ) banner();

  // Fix settings to be used for subobjects.
  settingsPtr->mode("Next:numberCount", 0);
  settingsPtr->mode("Next:numberShowLHA", 0);
  settingsPtr->mode("Next:numberShowInfo", 0);
  settingsPtr->mode("Next:numberShowProcess", 0);
  settingsPtr->mode("Next:numberShowEvent", 0);
  settingsPtr->flag("HadronLevel:all", false);
  vetoPrimaryProcess = clearSoftQCDFlags(*settingsPtr);
  settingsPtr->wvec("Init:plugins", {});

  // Create Pythia subobjects.
  for ( int i = MBIAS; i < ALL; ++i ) {
    pythia[i] = new Pythia(*settingsPtr, *particleDataPtr, false);
    pythia[i]->settings.mode("HeavyIon:mode", 1);
    pythia[i]->settings.flag("Beams:allowVertexSpread", false);
    if (i != MBIAS)
      pythia[i]->settings.mode("MultipartonInteractions:reuseInit", 0);
  }

  // Allow for user to override with a custom HIUserHooks.
  if ( HIHooksPtr ) HIHooksPtr->init(idProj, idTarg);

  // Initialize kinematics and cross sections.
  setBeamKinematics(idProjP, idTargP);
  for ( int i = MBIAS; i < ALL; ++i ) {
    if ( !pythia[i] ) continue;
    pythia[i]->settings.mode("Beams:frameType", 1);
    pythia[i]->settings.parm("Beams:eCM", beamSetupPtr->eCM);
  }

  // Initialize subobject for minimum bias processes.
  clearProcessLevel(*pythia[MBIAS]);
  pythia[MBIAS]->settings.flag("SoftQCD:all", true);
  pythia[MBIAS]->settings.mode("Beams:idA", idProjP);
  pythia[MBIAS]->settings.mode("Beams:idB", idTargP);
  if ( beamSetupPtr->frameType > 3 ) {
    pythia[MBIAS]->settings.parm("Beams:eA", beamSetupPtr->eA);
    pythia[MBIAS]->settings.parm("Beams:eB", beamSetupPtr->eB);
    pythia[MBIAS]->settings.mode("Beams:frameType", 2);
  }
  doLowEnergy = flag("LowEnergyQCD:all");
  if ( doLowEnergy )
    pythia[MBIAS]->settings.flag("LowEnergyQCD:all", true);

  pythia[MBIAS]->addUserHooksPtr(selectMB);
  init(MBIAS, "minimum bias processes");

  settingsPtr->wvec("Init:reuseMPIiDiffSys0",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys0"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys1",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys1"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys2",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys2"));
  settingsPtr->wvec("Init:reuseMPIiDiffSys3",
                    pythia[MBIAS]->settings.wvec("Init:reuseMPIiDiffSys3"));

  // Initialize semi-inclusive cross sections.
  sigTotNN.init();
  if ( !doLowEnergy && !sigTotNN.calc(idProjP, idTargP, beamSetupPtr->eCM))
    return false;

  // Set up nucleus geometry (projectile).
  if (HIHooksPtr && HIHooksPtr->hasProjectileModel())
    projPtr = HIHooksPtr->projectileModel();
  else
    projPtr = NucleusModel::create(mode("Angantyr:NucleusModelA"));
  if (!projPtr) {
    loggerPtr->ABORT_MSG("nucleus model not found for projectile");
    return false;
  }
  projPtr->initPtr(idProj, true, *infoPtr);
  if (!projPtr->init()) {
    loggerPtr->ABORT_MSG("projectile nucleus model failed to initialize");
    return false;
  }
  projPtr->setPN(beamSetupPtr->pAinit);

  // Set up nucleus geometry (target).
  if (HIHooksPtr && HIHooksPtr->hasTargetModel())
    targPtr = HIHooksPtr->targetModel();
  else
    targPtr = NucleusModel::create(mode("Angantyr:NucleusModelB"));
  if (!targPtr) {
    loggerPtr->ABORT_MSG("nucleus model not found for target");
    return false;
  }
  targPtr->initPtr(idTarg, false, *infoPtr);
  if (!targPtr->init()) {
    loggerPtr->ABORT_MSG("target nucleus model failed to initialize");
    return false;
  }
  targPtr->setPN(beamSetupPtr->pBinit);

  // Set up subcollision model.
  if ( HIHooksPtr && HIHooksPtr->hasSubCollisionModel() )
    collPtr = HIHooksPtr->subCollisionModel();
  else
    collPtr = SubCollisionModel::create(mode("Angantyr:CollisionModel"));
  if (!collPtr) {
    loggerPtr->ABORT_MSG("subcollision model not found");
    return false;
  }
  collPtr->initPtr(*projPtr, *targPtr, *pythia[MBIAS]->info.sigmaTotPtr,
                   *settingsPtr, *infoPtr, *rndmPtr);

  if ( doLowEnergy ) {
    lowEnergyCollPtr = SubCollisionModel::create(0);
    lowEnergyCollPtr->initPtr(*projPtr, *targPtr,
                              *pythia[MBIAS]->info.sigmaTotPtr,
                              *settingsPtr, *infoPtr, *rndmPtr);
    lowEnergyCollPtr->init(idProjP, idTargP, beamSetupPtr->eCM);
    lowEnergyCollPtr->initLowEnergy(pythia[MBIAS]->info.sigmaCmbPtr);
  }

  if (!collPtr->init(idProjP, idTargP, beamSetupPtr->eCM)) {
    loggerPtr->ABORT_MSG("subcollision model failed to initialize");
    return false;
  }
  if ( doLowEnergy && beamSetupPtr->eCM < eCMlow )
    hiInfo.avNDbSave = lowEnergyCollPtr->avNDB();
  else hiInfo.avNDbSave = collPtr->avNDB();

  // Set up impact parameter generator.
  if ( HIHooksPtr && HIHooksPtr->hasImpactParameterGenerator() )
    bGenPtr = HIHooksPtr->impactParameterGenerator();
  else
    bGenPtr = make_shared<ImpactParameterGenerator>();
  bGenPtr->initPtr(*infoPtr, *collPtr, *projPtr, *targPtr);
  if ( !bGenPtr->init() ) {
    loggerPtr->ABORT_MSG("impact parameter generator failed to initialize");
    return false;
  }

  if ( sabsMode >= 0 ) {
    // Initialize subobject for secondary absorptive processes.
    clearProcessLevel(*pythia[SASD]);
    Settings & sdabsopts = pythia[SASD]->settings;
    sdabsopts.flag("SoftQCD:singleDiffractive", true);

    setupSpecials(*pythia[SASD], "HIDiffraction:");
    setupSpecials(*pythia[SASD], "HIMultipartonInteractions:");
    setupSpecials(*pythia[SASD], "HIPDF:");
    setupSpecials(*pythia[SASD], "HISigmaDiffractive:");
    setupSpecials(*pythia[SASD], "HIBeamRemnants:");
    if ( sdabsopts.mode("Angantyr:SASDmode") > 0 ) {
      double pT0Ref = sdabsopts.parm("MultipartonInteractions:pT0Ref");
      double ecmRef = sdabsopts.parm("MultipartonInteractions:ecmRef");
      double ecmPow = sdabsopts.parm("MultipartonInteractions:ecmPow");
      double ecm = beamSetupPtr->eCM;
      sdabsopts.parm("Beams:eCM", ecm);
      double pT0     = pT0Ref * pow(ecm / ecmRef, ecmPow);
      sdabsopts.parm("MultipartonInteractions:pT0Ref", pT0, true);
      sdabsopts.parm("MultipartonInteractions:ecmRef", ecm, true);
      sdabsopts.parm("MultipartonInteractions:ecmPow", 0.0, true);
      sdabsopts.word("PDF:PomSet", "11");
      int reuseMpi = settingsPtr->mode("HeavyIon:SasdMpiReuseInit");
      if (reuseMpi != 0) {
        string initFile = settingsPtr->word("HeavyIon:SasdMpiInitFile");
        sdabsopts.mode("MultipartonInteractions:reuseInit", reuseMpi);
        sdabsopts.word("MultipartonInteractions:initFile", initFile);
        sdabsopts.wvec("Init:reuseMPIiDiffSys0",
                       settingsPtr->wvec("Init:reuseSasdMPIiDiffSys0"));
        sdabsopts.wvec("Init:reuseMPIiDiffSys1",
                       settingsPtr->wvec("Init:reuseSasdMPIiDiffSys1"));
        sdabsopts.wvec("Init:reuseMPIiDiffSys2",
                       settingsPtr->wvec("Init:reuseSasdMPIiDiffSys2"));
        sdabsopts.wvec("Init:reuseMPIiDiffSys3",
                     settingsPtr->wvec("Init:reuseSasdMPIiDiffSys3"));
      }
      if ( sdabsopts.mode("Angantyr:SASDmode") == 2 ) {
        sdabsopts.parm("Diffraction:mRefPomP", ecm);
        double sigND = sigTotNN.sigmaND();
        double mmin = sdabsopts.parm("Diffraction:mMinPert");
        double powp = sdabsopts.parm("HIDiffraction:mPowPomP");
        sdabsopts.parm("Diffraction:mPowPomP", powp, true);
        if ( powp > 0.0 ) sigND /= ((1.0 - pow(mmin/ecm, powp))/powp);
        else sigND /= log(ecm/mmin);
        sdabsopts.parm("Diffraction:sigmaRefPomP", sigND, true);
      }
      if ( sdabsopts.mode("Angantyr:SASDmode") >= 3 ) {
        sdabsopts.parm("Diffraction:mRefPomP", ecm);
        double sigND = sigTotNN.sigmaND();
        sdabsopts.parm("Diffraction:sigmaRefPomP", sigND, true);
        sdabsopts.parm("Diffraction:mPowPomP", 0.0);
      }
    }
    sdabsopts.mode("Beams:idA", idProjP);
    sdabsopts.mode("Beams:idB", idTargP);
    if ( beamSetupPtr->frameType > 3 ) {
      sdabsopts.parm("Beams:eA", beamSetupPtr->eA);
      sdabsopts.parm("Beams:eB", beamSetupPtr->eB);
      sdabsopts.mode("Beams:frameType", 2);
    }

    pythia[SASD]->addUserHooksPtr(selectSASD);

    init(SASD, "secondary absorptive processes as single diffraction.");

    settingsPtr->wvec("Init:reuseSasdMPIiDiffSys0",
                      sdabsopts.wvec("Init:reuseMPIiDiffSys0"));
    settingsPtr->wvec("Init:reuseSasdMPIiDiffSys1",
                      sdabsopts.wvec("Init:reuseMPIiDiffSys1"));
    settingsPtr->wvec("Init:reuseSasdMPIiDiffSys2",
                      sdabsopts.wvec("Init:reuseMPIiDiffSys2"));
    settingsPtr->wvec("Init:reuseSasdMPIiDiffSys3",
                    sdabsopts.wvec("Init:reuseMPIiDiffSys3"));
  }

  // Initialize subobject for hadronization.
  clearProcessLevel(*pythia[HADRON]);
  pythia[HADRON]->settings.flag("ProcessLevel:all", false);
  pythia[HADRON]->settings.flag("PartonLevel:all", false);
  pythia[HADRON]->settings.flag("HadronLevel:all", doHadronLevel);
  pythia[HADRON]->settings.mode("Beams:idA", idProj);
  pythia[HADRON]->settings.mode("Beams:idB", idTarg);
  pythia[HADRON]->settings.flag("LowEnergyQCD:all", false);

  // Initialize subobjects for signal processes.
  pythia[SIGPP]->settings.mode("Beams:idA", idProjP);
  pythia[SIGPP]->settings.mode("Beams:idB", idTargP);
  if ( idTargN ) {
    pythia[SIGPN]->settings.mode("Beams:idA", idProjP);
    pythia[SIGPN]->settings.mode("Beams:idB", idTargN);
  }
  if ( idProjN ) {
    pythia[SIGNP]->settings.mode("Beams:idA", idProjN);
    pythia[SIGNP]->settings.mode("Beams:idB", idTargP);
  }
  if ( idProjN && idTargN ) {
    pythia[SIGNN]->settings.mode("Beams:idA", idProjN);
    pythia[SIGNN]->settings.mode("Beams:idB", idTargN);
  }

  hasSignal = pythia[SIGPP]->settings.hasHardProc() ||
      pythia[SIGPP]->settings.mode("Beams:frameType") >= 4;
  if ( hasSignal ) {
    init(SIGPP, "signal process (pp)", 10);
    if ( idTargN ) init(SIGPN, "signal process (pn)", 10);
    if ( idProjN ) init(SIGNP, "signal process (np)", 10);
    if ( idProjN && idTargN ) init(SIGNN, "signal process (nn)", 10);
  }

  if (doHadronLevel) {
    if ( print )
      cout << " Angantyr Info: Initializing hadronisation processes." << endl;
  }
  settingsPtr->flag("ProcessLevel:all", false);
  infoPtr->hiInfo = &hiInfo;

  if ( cascadeMode ) updateMedium();

  return true;

}

//--------------------------------------------------------------------------

// Initialize a specific Pythia object and optionally run a number
// of events to get a handle of the cross section.

bool Angantyr::init(PythiaObject sel, string name, int n) {
  bool print = flag("HeavyIon:showInit") && !flag("Print:quiet");
  shared_ptr<InfoGrabber> ihg = make_shared<InfoGrabber>();
  pythia[sel]->addUserHooksPtr(ihg);
  if ( print ) cout << " Angantyr Info: Initializing " << name << "." << endl;
  if ( !pythia[sel]->init() ) return false;
  info[sel] = ihg->getInfo();
  if ( n <= 0 ) return true;
  if ( print ) cout << "Generating a few signal events for " << name
                    << " to build up statistics" << endl;
  for ( int i = 0; i < n; ++i ) pythia[sel]->next();
  return true;
}


//--------------------------------------------------------------------------

// Generate events and return EventInfo objects for different process
// types.

EventInfo Angantyr::getSignal(const SubCollision & coll) {
  if ( !hasSignal ) return EventInfo();
  int pytsel = SIGPP + coll.nucleons();
  int itry = MAXTRY;
  while ( itry-- ) {
    if ( pythia[pytsel]->next() ) {
      if ( pythia[pytsel]->event[0].pAbs2() != 0.0 )
        pythia[pytsel]->event.rotbst(toCMframe(pythia[pytsel]->event[1].p(),
                                               pythia[pytsel]->event[2].p()));
      return mkEventInfo(*pythia[pytsel], *info[pytsel], &coll);
    }
  }
  loggerPtr->WARNING_MSG("could not setup signal sub-collision");
  return EventInfo();
}

EventInfo Angantyr::getMBIAS(const SubCollision * coll, int procid) {
  int itry = MAXTRY;
  double bp = bMode && procid == 101? coll->bp: -1.0;
  double olapp = bMode < 0 && procid == 101? coll->olapp: -1.0;
  HoldProcess hold(selectMB, procid, bp, olapp);
  if ( allowIDAswitch )
    pythia[MBIAS]->setBeamIDs(coll->proj->id(), coll->targ->id());
  while ( --itry ) {
    if ( !pythia[MBIAS]->next(procid%10) ) continue;
    if (pythia[MBIAS]->info.code()%10 != procid%10) {
      loggerPtr->ERROR_MSG("MBIAS info code not equal to set procid",
                          "contact the authors");
      doAbort = true;
    }
    return mkEventInfo(*pythia[MBIAS], *info[MBIAS], coll);
  }
  return EventInfo();
}

EventInfo Angantyr::getSASD(const SubCollision * coll,
                            int procid, EventInfo * evp) {
  if ( sabsMode < 0 ) return getSABS(coll, procid, evp);
  int itry = MAXTRY;
  double bp = abs(bMode) > 1? coll->bp: -1.0;
  double olapp = bMode < -1? coll->olapp: -1.0;
  HoldProcess hold(selectSASD, procid, bp, olapp);
  if ( allowIDAswitch )
    pythia[MBIAS]->setBeamIDs(coll->proj->id(), coll->targ->id());
  while ( itry-- ) {
    if ( !pythia[SASD]->next(procid%10) ) continue;
    if (pythia[SASD]->info.code()%10 != procid%10) {
      loggerPtr->ERROR_MSG("SASD info code not equal to set procid",
                          "contact the authors");
      doAbort = true;
    }
    return mkEventInfo(*pythia[SASD], *info[SASD], coll);
  }
  return EventInfo();
}

EventInfo Angantyr::getSABS(const SubCollision * coll,
                            int procidIn, EventInfo * evp) {
  // Which side is to be diffractively excited.
  bool projex = ( procidIn == 103 );
  Pythia & mbgen = *pythia[MBIAS];
  // Get info about the sub collision
  int itry = MAXTRY;
  double eCM = beamSetupPtr->eCM;
  double bp = abs(bMode) > 1? coll->bp: -1.0;
  double olapp = bMode < -1? coll->olapp: -1.0;
  double x0 = pow2(sabsMinMX/eCM);
  double xc = pow2(sabsCutMX/eCM);

  // Get limits on the energy available for usage in the main event.
  double xlim = 1.0;
  if ( sabsMode < -2 ) {
    Vec4 prem = {};
    for ( int i : (projex? evp->targRems[coll->targ]:
                           evp->projRems[coll->proj]) )
      prem += evp->event[i].p();
    xlim = (projex? prem.pNeg()/eCM: prem.pPos()/eCM);
  }
  double e = -sabsEps;
  double xlime = pow(xlim, e);
  double x0e = pow(x0, e);

  if ( allowIDAswitch )
    mbgen.setBeamIDs(coll->proj->id(), coll->targ->id());
  double xmax = 0.0;
  int procid = procidIn;
  bool resample = true;
  while ( --itry ) {
    if ( resample ) {
      // Draw xmax from a 1/x^e distribution with special procedure for e==0.
      double R = rndmPtr->flat();
      xmax = ( e == 0.0 ? xlim*pow(x0/xlim, R):
               pow(xlime + R*(x0e - xlime), 1.0/e) );
      resample = false;
    }
    // Only do the special non-diffractive event if xmax > xc,
    // otherwise we do a standard single diffractive event
    // (specified in procidIn).
    procid = procidIn;
    if ( xmax > xc ) procid = 101;
    else xmax = xc;
    if ( projex ) xmax = -xmax;

    // Setup the selectMB UserHook temporarily (is reset again when
    // "hold" is out of scope).
    HoldProcess hold(selectMB, procid, bp, olapp, xmax);
    if ( !mbgen.next(procid % 10) ) {
      if ( selectMB->failcount >= 100 ) resample = true;
      continue;
    }
    // LowEnergyQCD (151) processes cannot be handled here. so try
    // again until we get a SoftQCD one (101).
    if ( mbgen.info.code() == 151 ) continue;
    if (mbgen.info.code()%10 != procid%10) {
      loggerPtr->ERROR_MSG("SABS info code not equal to set procid",
                           "contact the authors");
      doAbort = true;
    }
    if ( procid%10 == 1 &&
         !fixSecondaryAbsorptive(mbgen.event, xmax) )
      continue;
    return mkEventInfo(mbgen, *info[MBIAS], coll);
  }
  return EventInfo();
}

//--------------------------------------------------------------------------

// Generate primary absorptive (non-diffractive) nucleon-nucleon
// sub-collisions.

bool Angantyr::genAbs(SubCollisionSet& subCollsIn,
  list<EventInfo>& subEventsIn) {
  // The fully absorptive.
  vector<const SubCollision*> abscoll;
   // The partly absorptive.
  vector<const SubCollision*> abspart;
  // The non-diffractive and signal events.
  multiset<EventInfo> ndeve, sigeve;

  // Select the primary absorptive sub collisions.
  for (const SubCollision& subColl : subCollsIn) {

    if ( subColl.type != SubCollision::ABS ) continue;
    if (!subColl.proj->done() && !subColl.targ->done() ) {
      abscoll.push_back(&subColl);
      if ( bMode != 0 ) {
        EventInfo ie = getMBIAS(&subColl, 101);
        if (ie.code%10  != 1) {
          loggerPtr->ERROR_MSG("ND code not equal to 101",
                            "contact the authors");
          doAbort = true;
        }
        ndeve.insert(ie);
      }
      subColl.proj->select();
      subColl.targ->select();
    } else
      abspart.push_back(&subColl);
  }

  if ( abscoll.empty() ) return true;

  int Nabs = abscoll.size();
  int Nadd = abspart.size();

  if ( bMode == 0 ) {
    for ( int i = 0; i < Nabs + Nadd; ++i ) {
      EventInfo ie = getMBIAS(i < Nabs? abscoll[i]: abspart[i - Nabs], 101);
      if (ie.code%10 != 1) {
        loggerPtr->ERROR_MSG("ND code not equal to 101",
                            "contact the authors");
        doAbort = true;
      }
      ndeve.insert(ie);
    }
  }
  vector<int> Nii(4, 0);
  vector<double> w(4, 0.0);
  double wsum = 0.0;
  double P1 = 1.0;
  if ( hasSignal ) {

    // Count how many potential absorpitve collisions there are for
    // each iso-spin combination.
    for ( int i = 0, N = abscoll.size(); i < N; ++i )
      ++Nii[abscoll[i]->nucleons()];
    for ( int i = 0, N = abspart.size(); i < N; ++i )
      ++Nii[abspart[i]->nucleons()];

    if ( Nii[0] )
      w[0] = pythia[SIGPP]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[1] )
      w[1] = pythia[SIGPN]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[2] )
      w[2] = pythia[SIGNP]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();
    if ( Nii[3] )
      w[3] = pythia[SIGNN]->info.sigmaGen()*MB2FMSQ/collPtr->sigND();

    wsum = Nii[0]*w[0] + Nii[1]*w[1] + Nii[2]*w[2] + Nii[3]*w[3];
    P1 = 1.0 - pow(1.0 - w[0], Nii[0])*pow(1.0 - w[1], Nii[1])*
               pow(1.0 - w[2], Nii[2])*pow(1.0 - w[3], Nii[3]);

  }

  bool noSignal = hasSignal;

  // *** THINK *** Is it ok to always pair the hardest events with the
  // *** most central sub-collisions, or will this introduce a strange
  // *** bias?
  multiset<EventInfo>::iterator it = ndeve.begin();
  EventInfo ei;
  for ( int i = 0, N = abscoll.size(); i < N; ++i ) {
    int b = abscoll[i]->nucleons();
    if ( Nii[b]
         && ( noSignal || w[b]*(wsum/P1 - 1.0)/(wsum - w[b]) > rndmPtr->flat())
         && (ei = getSignal(*abscoll[i])).ok ) {
      noSignal = false;
    }
    else
      ei =*it++;
    subEventsIn.push_back(ei);
    if ( !setupFullCollision(subEventsIn.back(), *abscoll[i],
                             Nucleon::ABS, Nucleon::ABS) )
      return false;
  }

  if ( noSignal ) return false;

  hiInfo.reweight(P1);

  return true;

}

//--------------------------------------------------------------------------

// Add secondary absorptive sub-collisions to the primary ones.

void Angantyr::addSASD(const SubCollisionSet& subCollsIn) {
  // Collect absorptively wounded nucleons in secondary
  // sub-collisions.
  int ntry = mode("Angantyr:SDTries");
  if ( settingsPtr->isMode("HI:SDTries") )
    ntry = mode("HI:SDTries");
  for (const SubCollision& subColl : subCollsIn)
    if ( subColl.type == SubCollision::ABS ) {
      if ( subColl.targ->done() && !subColl.proj->done() ) {
        EventInfo * evp = subColl.targ->event();
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getSASD(&subColl, 103, evp);
          if ( addNucleonExcitation(*evp, add, true) ) {
            subColl.proj->select(*evp, Nucleon::ABS);
            break;
          }
          if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
        }
      } else if ( subColl.proj->done() && !subColl.targ->done() ) {
        EventInfo * evp = subColl.proj->event();
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getSASD(&subColl, 104, evp);
          if ( addNucleonExcitation(*evp, add, true) ) {
            subColl.targ->select(*evp, Nucleon::ABS);
            break;
          }
          if ( itry == ntry - 1 ) hiInfo.failedExcitation(subColl);
        }
      }
    }
}

//--------------------------------------------------------------------------

// Add primary sub-collisions, of the given types.

bool Angantyr::addSubCollisions(const SubCollisionSet& subCollsIn,
  list<EventInfo> & subEventsIn, vector<CollDesc> colldescs) {
  using std::get;
  // Collect full single diffraction collisions.
  for (const SubCollision& subColl : subCollsIn)
    if ( !subColl.proj->done() && !subColl.targ->done() ) {
      for ( auto desc : colldescs ) {
        if ( subColl.type == get<0>(desc) ) {
          subEventsIn.push_back(getMBIAS(&subColl, get<1>(desc)));
          if ( !setupFullCollision(subEventsIn.back(), subColl,
                                   get<2>(desc), get<3>(desc)) )
            return false;
        }
      }
    }

  return true;
}

//--------------------------------------------------------------------------

// Add secondary sub-collisions of given types to primary ones.

void Angantyr::addSecondaries(const SubCollisionSet& subCollsIn,
  vector<CollDesc> colldescs) {
  using std::get;
  int ntry = mode("Angantyr:SDTries");
  if ( settingsPtr->isMode("HI:SDTries") )  ntry = mode("HI:SDTries");
  for (const SubCollision& subColl : subCollsIn) {
    bool isExc = ( subColl.type == SubCollision::DDE ||
                   subColl.type == SubCollision::SDEP ||
                   subColl.type == SubCollision::SDET ||
                   subColl.type == SubCollision::LEXC );
    for ( auto desc : colldescs ) {
      if ( !subColl.proj->done() && subColl.type == get<0>(desc) &&
           get<2>(desc) != Nucleon::UNWOUNDED ) {
        EventInfo * evp = subColl.targ->event();
        if ( evp->code == 158 || evp->code == 159 ) continue;
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getMBIAS(&subColl, get<1>(desc));
          if ( addNucleonExcitation(*evp, add, false) ) {
            subColl.proj->select(*evp, get<2>(desc));
            break;
          }
          if ( itry == ntry - 1 && isExc ) hiInfo.failedExcitation(subColl);
        }
      }
      if ( !subColl.targ->done() && subColl.type == get<0>(desc) &&
           get<3>(desc) != Nucleon::UNWOUNDED ) {
        EventInfo * evp = subColl.proj->event();
        if ( evp->code == 158 || evp->code == 159 ) continue;
        for ( int itry = 0; itry < ntry; ++itry ) {
          EventInfo add = getMBIAS(&subColl, get<1>(desc));
          if ( addNucleonExcitation(*evp, add, false) ) {
            subColl.targ->select(*evp, get<3>(desc));
            break;
          }
          if ( itry == ntry - 1 && isExc ) hiInfo.failedExcitation(subColl);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------

// Shift an event in impact parameter from the nucleon-nucleon
// sub-collision to the overall nucleus-nucleus frame. It is assumed
// that all partonic vertices are given in units of femtometers.

EventInfo & Angantyr::shiftEvent(EventInfo & ei) {
  if ( HIHooksPtr && HIHooksPtr->canShiftEvent() )
    return HIHooksPtr->shiftEvent(ei);

  double ymax = ei.event[1].y();
  Vec4 bmax = ei.coll->proj->bPos();
  double ymin = ei.event[2].y();
  Vec4 bmin = ei.coll->targ->bPos();
  for ( int i = 0, N = ei.event.size(); i < N; ++i ) {
    Vec4 shift = bmin + (bmax - bmin)*(ei.event[i].y() - ymin)/(ymax - ymin);
    ei.event[i].vProdAdd( shift * FM2MM);
  }
  return ei;
}

//--------------------------------------------------------------------------

// Prepare a primary sub-collision.

bool Angantyr::
setupFullCollision(EventInfo & ei, const SubCollision & coll,
                   Nucleon::Status projStatus, Nucleon::Status targStatus) {
  if ( !ei.ok ) return false;
  coll.proj->select(ei, projStatus);
  coll.targ->select(ei, targStatus);
  ei.coll = &coll;
  ei.projs.clear();
  ei.projs[coll.proj] = make_pair(1, ei.event.size());
  ei.targs.clear();
  ei.targs[coll.targ] = make_pair(2, ei.event.size());
  shiftEvent(ei);
  // Note that -203 is a special status code for identifying incoming
  // particles in a sub-event.
  ei.event[1].status(-203);
  ei.event[1].mother1(1);
  ei.event[1].mother2(0);
  ei.event[2].status(-203);
  ei.event[2].mother1(2);
  ei.event[2].mother2(0);
  return fixIsoSpin(ei);
}

//--------------------------------------------------------------------------

// Trace a particle back to one of the beams in an event.

int Angantyr::getBeam(Event & ev, int i) {
  if ( int mom = ev[i].mother1() ) {
    if ( ev[mom].status() != -203 && ev[mom].mother1() < mom )
      return getBeam(ev, mom);
    else
      return mom;
  }
  else
    return i;
}

//--------------------------------------------------------------------------

// Minimum-bias sub-collisions are always generated as p-p events.
// We assume the generator is isospin invariant, and
// post–hoc flip the beam IDs and remnant quark content to emulate
// p-n, n-p, and n-n collisions.

bool Angantyr::fixIsoSpin(EventInfo & ei) {
  if ( HIHooksPtr && HIHooksPtr->canFixIsoSpin() )
    return HIHooksPtr->fixIsoSpin(ei);

  // Check if isospin needs fixing.
  int pshift = 0, tshift = 0;
  if ( ei.event[1].id() == 2212 && ei.coll->proj->id() == 2112 )
    pshift = 1;
  if ( ei.event[1].id() == -2212 && ei.coll->proj->id() == -2112 )
    pshift = -1;
  if ( pshift )
    ei.event[1].id(pshift*2112);
  if ( ei.event[2].id() == 2212 && ei.coll->targ->id() == 2112 )
    tshift = 1;
  if ( ei.event[2].id() == -2212 && ei.coll->targ->id() == -2112 )
    tshift = -1;
  if ( tshift )
    ei.event[2].id(tshift*2112);

  if ( !pshift && !tshift ) return true;

  // Try to find corresponding remnants that change flavour
  for ( int i = ei.event.size()  - 1; i > 2 && ( pshift || tshift ); --i ) {
    if ( pshift && ( isRemnant(ei, i) || ei.event[i].status() == 14 )
         &&  getBeam(ei.event, i) == 1 ) {
      int newid = 0;
      if ( ei.event[i].id() == 2*pshift ) newid = 1*pshift;
      if ( ei.event[i].id() == 2101*pshift ) newid = 1103*pshift;
      if ( ei.event[i].id() == 2103*pshift ) newid = 1103*pshift;
      if ( ei.event[i].id() == 2203*pshift ) newid = 2103*pshift;
      if ( ei.event[i].id() == 2212*pshift ) newid = 2112*pshift;
      if ( newid ) {
        ei.event[i].id(newid);
        pshift = 0;
        continue;
      }
    }
    if ( tshift && ( isRemnant(ei, i) || ei.event[i].status() == 14 )
         &&  getBeam(ei.event, i) == 2 ) {
      int newid = 0;
      if ( ei.event[i].id() ==    2*tshift ) newid =    1*tshift;
      if ( ei.event[i].id() == 2101*tshift ) newid = 1103*tshift;
      if ( ei.event[i].id() == 2103*tshift ) newid = 1103*tshift;
      if ( ei.event[i].id() == 2203*tshift ) newid = 2103*tshift;
      if ( ei.event[i].id() == 2212*tshift ) newid = 2112*tshift;
      if ( newid ) {
        ei.event[i].id(newid);
        tshift = 0;
        continue;
      }
    }
  }

  if ( !pshift && !tshift ) return true;

  // Try to find any final state quark that we modify, preferably far
  // in the beam direction.
  // For projectile (beam 1), pick the most forward u-quark (largest y).
  // For target (beam 2), pick the most backward u-quark (smallest y).
  int qselp = 0;
  int qselt = 0;
  double yselp = 0.0;
  double yselt = 0.0;
  for ( int i = ei.event.size()  - 1; i > 2 && ( pshift || tshift ); --i ) {
    if ( pshift && ei.event[i].isFinal() && ei.event[i].id() == 2*pshift) {
      if ( ei.event[i].y() > yselp ) {
        qselp = i;
        yselp = ei.event[i].y();
      }
    }
    if ( tshift && ei.event[i].isFinal() && ei.event[i].id() == 2*tshift) {
      if ( ei.event[i].y() < yselt ) {
        qselt = i;
        yselt = ei.event[i].y();
      }
    }
  }
  if ( qselp ) {
    ei.event[qselp].id(1*pshift);
    pshift = 0;
  }
  if ( qselt ) {
    ei.event[qselt].id(1*tshift);
    tshift = 0;
  }

  return !pshift && !tshift;

}

//--------------------------------------------------------------------------

// Find recoilers in a primary sub-collisions to conserve energy and
// momentum when adding a secondary one. Not actually used yet.

vector<int> Angantyr::
findRecoilers(const Event & e, bool tside, int beam, int end,
              const Vec4 & pdiff, const Vec4 & pbeam) {
  vector<int> ret;
  multimap<double,int> ordered;
  double mtd2 = pdiff.m2Calc() + pdiff.pT2();
  int dir = tside? -1: 1;
  double ymax = -log(pdiff.pNeg());
  if ( tside ) ymax = -log(pdiff.pPos());
  for ( int i = beam, N = end; i < N; ++i )
    if ( e[i].status() > 0 )
      ordered.insert(make_pair(e[i].y()*dir, i));
  Vec4 prec;
  double pzfree2 = 0.0;
  multimap<double,int>::iterator it = ordered.begin();
  while ( it != ordered.end() ) {
    if ( it->first > ymax ) break;
    int i = (*it++).second;
    Vec4 test = prec + e[i].p();
    double mtr2 = test.m2Calc() + test.pT2();
    double S = (pbeam + test).m2Calc();
    double pz2 = 0.25*(pow2(S - mtr2 - mtd2) - 4.0*mtr2*mtd2)/S;
    if ( pz2 < pzfree2 ) break;
    prec = test;
    pzfree2 = pz2;
    ret.push_back(i);
  }

  // *** THINK! *** Is this the best way?
  return ret;

}

//--------------------------------------------------------------------------

// Add a secondary sub-collision to a primary one.

bool Angantyr::addNucleonExcitation(EventInfo & ei, EventInfo & sub,
                                    bool colConnect) {
  if ( !sub.ok ) return false;
  if ( sabsMode < 0 || mode("Angantyr:testMode") > 0 )
    return addNucleonExcitation2(ei, sub, colConnect);
  if ( mode("Angantyr:testMode") < -1 )
    cout << "Adding " << sub.code << " to " << ei.code << endl;
  fixIsoSpin(sub);
  shiftEvent(sub);

  if ( HIHooksPtr && HIHooksPtr->canAddNucleonExcitation() )
    return HIHooksPtr->addNucleonExcitation(ei, sub, colConnect);

  typedef map<Nucleon *, pair<int, int> >::iterator NucPos;
  bool tside = false;
  NucPos recnuc = ei.projs.find(sub.coll->proj);
  if ( recnuc != ei.projs.end() ) tside = true;
  NucPos rectarg = ei.targs.find(sub.coll->targ);
  if ( rectarg != ei.targs.end() ) {
    if ( tside ) loggerPtr->WARNING_MSG("nucleon already added");
    tside = false;
    recnuc = rectarg;
  }

  // First get the projectile system to take recoil and their momentum.
  int olddiff = tside? 4: 3;
  int beam = tside? 2: 1;
  int recbeam = recnuc->second.first;
  int recend = recnuc->second.second;
  Vec4 pbeam = sub.event[beam].p();
  Vec4 pdiff = sub.event[olddiff].p();
  if ( sub.code == 106 ) pdiff += sub.event[5].p();
  if ( sub.code == 101 )
    pdiff = sub.event[0].p() - sub.event[sub.event.size() - 1].p();
  if ( sub.code == 157 && tside )
    pdiff = sub.event[0].p() - sub.event[3].p();
  if ( sub.code == 157 && !tside )
    pdiff = sub.event[0].p() - sub.event[4].p();
  vector<int> rec;
  Vec4 prec;
  if ( HIHooksPtr && HIHooksPtr->canFindRecoilers() )
    rec = HIHooksPtr->findRecoilers(ei.event, tside, recbeam, recend,
                                    pdiff, pbeam);
  else if ( recoilerMode == 2 )
    rec = findRecoilers(ei.event, tside, recbeam, recend, pdiff, pbeam);
  else {
    if ( tside && ei.event[3].status() > 0 &&
         ( ei.code == 104 || ei.code == 157 ) )
      rec.push_back(3);
    else if ( !tside && ei.event[4].status() > 0 &&
      ( ei.code == 104 || ei.code == 157 ) )
      rec.push_back(3);
    else if ( tside && ei.event[3].status() > 0 &&
              ( ei.code == 102 || ei.code == 106 || ei.code == 157 ) )
      rec.push_back(3);
    else if ( !tside && ei.event[4].status() > 0 &&
              ( ei.code == 102 || ei.code == 106 || ei.code == 157 ) )
      rec.push_back(4);
    else
      for ( int i = recbeam, N = recend; i < N; ++i )
        if ( isRemnant(ei, i) && getBeam(ei.event, i) == recbeam )
          rec.push_back(i);
  }
  if ( rec.empty() ) return false;
  for ( int i = 0, N = rec.size(); i < N; ++i ) prec += ei.event[rec[i]].p();

  if ( mode("Angantyr:testMode") < -1 ) {
    for ( int i : rec ) cout << i << " ";
    cout << "-> " << pdiff;
  }

  // Find the transform to the recoilers and the diffractive combined cms.
  pair<RotBstMatrix,RotBstMatrix> R12;
  if ( !getTransforms(prec, pdiff, pbeam, R12) )
    return false;

  // Transform the recoilers.
  for ( int i = 0, N = rec.size(); i < N; ++i )
    ei.event[rec[i]].rotbst(R12.first);

  // Copy the event and transform and offset the particles
  // appropriately. Note that -203 is a special status code for
  // identifying incoming particles in a sub-event.
  int newbeam = ei.event.size();
  ei.event.append(sub.event[beam]);
  ei.event.back().status(-203);
  ei.event.back().mother1(beam);
  ei.event.back().mother2(0);
  ei.event.back().daughter1(ei.event.size());
  int newdiff = ei.event.size();
  int nextpos = 5;
  if ( sub.code == 101 || sub.code == 157 ) {
    ei.event.append(9902210, -15, 0, 0, pdiff, pdiff.mCalc());
    nextpos = 3;
    olddiff = tside? 1: 2;
  } else ei.event.append(sub.event[olddiff]);
  ei.event.back().rotbst(R12.second);
  ei.event.back().mother1(newbeam);
  ei.event.back().mother2(0);
  if ( sub.code == 102 ) {
    if ( tside )
      ei.targs[sub.coll->targ] = make_pair(newbeam, ei.event.size());
    else
      ei.projs[sub.coll->proj] = make_pair(newbeam, ei.event.size());
    if ( mode("Angantyr:testMode") < -1 ) cout << "worked" << endl;
    return true;
  }

  int idoff = tside? newdiff - olddiff: newdiff - olddiff - 1;
  if ( sub.code == 106 ) {
    // Special handling of central diffraction.
    ++newdiff;
    ++nextpos;
    idoff = newdiff - 5;
    ei.event.append(sub.event[5]);
    ei.event.back().rotbst(R12.second);
    ei.event.back().mother1(newbeam);
    ei.event.back().mother2(0);
  }
  ei.event.back().daughter1(sub.event[olddiff].daughter1() + idoff);
  if ( sub.code != 101 )
    ei.event.back().daughter2(sub.event[olddiff].daughter2() + idoff);

  int coloff = ei.event.lastColTag();
  // Add energy to zeroth line and calculate new invariant mass.
  ei.event[0].p( ei.event[0].p() + pbeam );
  ei.event[0].m( ei.event[0].mCalc() );
  for (int i = nextpos; i < sub.event.size(); ++i) {
    if ( sub.code == 101 && i + 1 == sub.event.size() ) break;
    Particle temp = sub.event[i];

    // Add offset to nonzero mother, daughter and colour indices.
    if ( temp.mother1() == olddiff ) temp.mother1(newdiff);
    else if ( temp.mother1() > 0 ) temp.mother1(temp.mother1() + idoff );
    if ( temp.mother2() == olddiff ) temp.mother2(newdiff);
    else if ( temp.mother2() > 0 ) temp.mother2(temp.mother2() + idoff );
    if ( temp.daughter1() > 0 ) temp.daughter1( temp.daughter1() + idoff );
    if ( temp.daughter2() > 0 ) temp.daughter2( temp.daughter2() + idoff );
    if ( temp.col() > 0 ) temp.col( temp.col() + coloff );
    if ( temp.acol() > 0 ) temp.acol( temp.acol() + coloff );
    temp.rotbst(R12.second);
    // Append particle to summed event.
    ei.event.append( temp );
  }

  addJunctions(ei.event, sub.event, coloff);

  if ( tside )
    ei.targs[sub.coll->targ] = make_pair(newbeam, ei.event.size());
  else
    ei.projs[sub.coll->proj] = make_pair(newbeam, ei.event.size());

  if ( mode("Angantyr:testMode") < -1 ) cout << "worked" << endl;
  return true;

}

// Add a secondary sub-collision to a primary one (new version).

bool Angantyr::addNucleonExcitation2(EventInfo & ei, EventInfo & sub,
                                     bool colConnect) {
  fixIsoSpin(sub);
  shiftEvent(sub);
  if ( mode("Angantyr:testMode") > 1 )
    cout << "Adding " << sub.code << " to " << ei.code << endl;
  if ( HIHooksPtr && HIHooksPtr->canAddNucleonExcitation() )
    return HIHooksPtr->addNucleonExcitation(ei, sub, colConnect);

  // Find out which nucleon in the new sub-event that is to be added.
  typedef map<Nucleon *, pair<int, int> >::iterator NucPos;
  bool tside = false;
  NucPos recnuc = ei.projs.find(sub.coll->proj);
  if ( recnuc != ei.projs.end() ) tside = true;
  NucPos rectarg = ei.targs.find(sub.coll->targ);
  if ( rectarg != ei.targs.end() ) {
    if ( tside ) {
      loggerPtr->WARNING_MSG("nucleon already added");
      return false;
    }
    tside = false;
    recnuc = rectarg;
  }

  // This is the beam to be added with its excited state.
  int beam = tside? 2: 1;
  Vec4 pbeam = sub.event[beam].p();
  int iqel = (tside? sub.projEl[sub.coll->proj]: sub.targEl[sub.coll->targ]);
  if ( !iqel ) return false;
  int rmbeam = 3 - beam;
  Vec4 ppom = sub.event[rmbeam].p() - sub.event[iqel].p();
  Vec4 pdiff = pbeam + ppom;

  // This is the beam for which the remnants will receive a shift in
  // momentum due to more energy is taken out to go into the new sub
  // event.
  int recbeam = recnuc->second.first;
  vector<int> & rec = ( tside? ei.projRems[sub.coll->proj]:
                               ei.targRems[sub.coll->targ] );
  if ( rec.empty() ) return false;
  Vec4 prec;
  for ( int i = 0, N = rec.size(); i < N; ++i ) prec += ei.event[rec[i]].p();

  if ( mode("Angantyr:testMode") > 1 ) {
    for ( int i : rec ) cout << i << " ";
    cout << "-> " << pdiff;
  }

  // Find the transform to the recoilers and the diffractive combined cms.
  pair<RotBstMatrix,RotBstMatrix> R12;
  if ( !getTransforms(prec, pdiff, pbeam, R12) ) return false;

  if ( mode("Angantyr:testMode") ) {
    // Just checking:
    Vec4 pbefore = prec + pbeam;
    Vec4 pafter = R12.first*prec + R12.second*pdiff;
    Vec4 pchange = pafter - pbefore;
    if ( pchange.pT() > 0.001 ) {
      cout << "There is some problem with the transforms!"
           << pbefore << pafter << pchange;
    }
  }

  // Copy and transform the recoilers.
  for ( int i = 0, N = rec.size(); i < N; ++i ) {
    rec[i] = ei.event.copy(rec[i], 66);
    ei.event[rec[i]].rotbst(R12.first);
  }

  // Now we start copying over the relevant parts of the new First
  // make room for the beams. Note that -203 is a special status code
  // for identifying incoming particles in a sub-event.
  int ioff = ei.event.size() - 1;
  auto shift = [&](int idx) ->  int { return idx > 0? idx + ioff: idx; };
  Particle newbeam = sub.event[beam];
  newbeam.status(-203);
  newbeam.mothers(beam);
  newbeam.daughters(shift(newbeam.daughter1()), shift(newbeam.daughter2()));
  Particle newpom(990, -13);
  newpom.mothers(recbeam);
  newpom.daughters(shift(sub.event[rmbeam].daughter1()),
                   shift(sub.event[rmbeam].daughter2()));
  newpom.p(R12.second*ppom);
  newpom.m(newpom.mCalc());
  ei.event.append(tside? newpom: newbeam);
  ei.event.append(tside? newbeam: newpom);

  // Add energy to zeroth line and calculate new invariant mass.
  ei.event[0].p( ei.event[0].p() + pbeam );
  ei.event[0].m( ei.event[0].mCalc() );

  // Copy everyting else to the main event.
  int coloff = ei.event.lastColTag();
  auto cshift = [&](int cidx) -> int { return cidx > 0? cidx + coloff: cidx; };
  for (int i = 3; i < sub.event.size(); ++i) {
    Particle p = sub.event[i];
    p.mothers(shift(p.mother1()), shift(p.mother2()));
    p.daughters(shift(p.daughter1()), shift(p.daughter2()));
    p.cols(cshift(p.col()), cshift(p.acol()));
    p.rotbst(R12.second);
    // Append particle to summed event.
    ei.event.append(p);
  }
  addJunctions(ei.event, sub.event, coloff);

  // Now remove the elastically scattered old beam, and setup info
  // about the new subevent in the main event.
  ei.event.remove(iqel + ioff, iqel + ioff);
  ei.event.remove(3 - beam + ioff, 3 - beam + ioff);
  vector<int> newrems = tside? sub.targRems[sub.coll->targ]
    : sub.projRems[sub.coll->proj];
  for ( int & i : newrems ) i = shift(i > iqel? i - 2: i - 1);
  int & newiqel = tside? sub.targEl[sub.coll->targ]
    : sub.projEl[sub.coll->proj];
  if ( newiqel ) newiqel = shift(newiqel > iqel? newiqel - 2: newiqel - 1);


  if ( tside ) {
    ei.targs[sub.coll->targ] = make_pair(1 + ioff, ei.event.size());
    ei.targRems[sub.coll->targ] = newrems;
    ei.targEl[sub.coll->targ] = newiqel;
  } else {
    ei.projs[sub.coll->proj] = make_pair(1 + ioff, ei.event.size());
    ei.projRems[sub.coll->proj] = newrems;
    ei.projEl[sub.coll->proj] = newiqel;
  }

  if ( mode("Angantyr:testMode") ) {
    Vec4 ptot;
    for ( auto p : ei.event ) if ( p.isFinal() ) ptot += p.p();
    if ( ptot.pT() > 0.001 )
      cout << "Energy imbalance when adding nuclon excitation:" << endl
           << ptot;
  }
  if ( mode("Angantyr:testMode") > 1 ) cout << "worked" << endl;

  return true;

}

//--------------------------------------------------------------------------

// Calculate boosts to shuffle momenta when adding secondary
// sub-collisions.

bool
Angantyr::getTransforms(Vec4 prec, Vec4 pdiff, const Vec4 & pbeam,
                      pair<RotBstMatrix,RotBstMatrix> & R12) {
  RotBstMatrix Ri;
  Ri.toCMframe(pbeam, prec);
  Vec4 pr1 = prec;
  Vec4 pb1 = pbeam;
  Vec4 pd1 = pdiff;
  pr1.rotbst(Ri);
  pb1.rotbst(Ri);
  pd1.rotbst(Ri);
  Vec4 pr2 = pr1;
  if ( pd1.pT() >= abs(pr2.pz()) ) {
    return false;
  }
  double the = asin(pd1.pT()/abs(pr2.pz()));
  RotBstMatrix R1;
  R1.rot(the, pd1.phi());
  pr2.rotbst(R1);

  double S = (prec + pbeam).m2Calc();
  double mtr2 = pr2.pT2() + pr2.m2Calc();
  double mtd2 = pd1.pT2() + pd1.m2Calc();
  if ( sqrt(S) <= sqrt(mtr2) + sqrt(mtd2) ) {
    return false;
  }
  double z2 = 0.25*(mtr2*mtr2 + (mtd2 - S)*(mtd2 - S) - 2.0*mtr2*(mtd2 + S))/S;
  if ( z2 <= 0.0 ) {
    return false;
  }
  double z = sqrt(z2);
  double ppo2 = pow2(pr2.pNeg());
  double ppn2 = pow2(z + sqrt(z2 + mtr2));
  R1.bst(0.0, 0.0, -(ppn2 - ppo2)/(ppn2 + ppo2));

  ppo2 = pow2(pd1.pPos());
  ppn2 = pow2(z + sqrt(z2 + mtd2));
  RotBstMatrix R2;
  R2.bst(0.0, 0.0, (ppn2 - ppo2)/(ppn2 + ppo2));
  Vec4 pr3 = pr1;
  pr3.rotbst(R1);
  Vec4 pd3 = pd1;
  pd3.rotbst(R2);

  RotBstMatrix Rf = Ri;
  Rf.invert();
  Vec4 pr4 = pr3;
  pr4.rotbst(Rf);
  Vec4 pd4 = pd3;
  pd4.rotbst(Rf);

  R12.first = R12.second = Ri;
  R12.first.rotbst(R1);
  R12.second.rotbst(R2);
  R12.first.rotbst(Rf);
  R12.second.rotbst(Rf);
  prec.rotbst(R12.first);
  pdiff.rotbst(R12.second);

  return true;

}

//--------------------------------------------------------------------------

// Add sub-events together taking special care with the status of the
// incoming nucleons, and also handle the junctions correctly.

void Angantyr::addSubEvent(Event & evnt, Event & sub) {

  int idoff = evnt.size() - 1;
  int coloff = evnt.lastColTag();

  for (int i = 1; i < sub.size(); ++i) {
    Particle temp = sub[i];

    // Add offset to nonzero mother, daughter and colour indices.
    if ( temp.status() == -203 ) {
      // Incoming particle to a sub-event
      temp.status(-13);
    } else {
      if ( temp.mother1() > 0 ) temp.mother1(temp.mother1() + idoff );
      if ( temp.mother2() > 0 ) temp.mother2( temp.mother2() + idoff );
    }
    if ( temp.daughter1() > 0 ) temp.daughter1( temp.daughter1() + idoff );
    if ( temp.daughter2() > 0 ) temp.daughter2( temp.daughter2() + idoff );
    if ( temp.col() > 0 ) temp.col( temp.col() + coloff );
    if ( temp.acol() > 0 ) temp.acol( temp.acol() + coloff );
    // Append particle to summed event.
    evnt.append( temp );
  }

  addJunctions(evnt, sub, coloff);

}

void Angantyr::addJunctions(Event & ev, Event & addev, int coloff) {

  // Read out junctions one by one.
  Junction tempJ;
  int begCol, endCol;
  for (int i = 0; i < addev.sizeJunction(); ++i) {
    tempJ = addev.getJunction(i);

    // Add colour offsets to all three legs.
    for (int  j = 0; j < 3; ++j) {
      begCol = tempJ.col(j);
      endCol = tempJ.endCol(j);
      if (begCol > 0) begCol += coloff;
      if (endCol > 0) endCol += coloff;
      tempJ.cols( j, begCol, endCol);
    }
    // Append junction to summed event.
    ev.appendJunction( tempJ );
  }
}

//--------------------------------------------------------------------------

// Special function to generate secondary absorptive events as single
// diffraction. Called from Angantyr::next() and used for debugging
// and tuning purposes.

bool Angantyr::nextSASD(int procid) {
  if ( sabsMode < 0 ) {
    loggerPtr->ERROR_MSG("SASD event requester but Angantyr:SASDmode < 0.");
    return false;
  }

  Nucleon dummy;
  double bp = pythia[SASD]->parm("Angantyr:SDTestB");
  SubCollision coll(dummy, dummy, bp*collPtr->avNDB(), bp, -1.0,
    SubCollision::ABS);
  EventInfo eidummy;
  EventInfo ei = getSASD(&coll, procid, &eidummy);
  if ( !ei.ok ) return false;
  pythia[HADRON]->event = ei.event;
  updateInfo();
  if (doHadronLevel) {
    if ( HIHooksPtr && HIHooksPtr->canForceHadronLevel() ) {
      if ( !HIHooksPtr->forceHadronLevel(*pythia[HADRON]) ) return false;
    } else {
      if ( !pythia[HADRON]->forceHadronLevel(false) ) return false;
    }
  }
  return true;
}

//--------------------------------------------------------------------------

// Reset the main event.

void Angantyr::resetEvent() {

  Event & etmp = pythia[HADRON]->event;
  unifyFrames();
  etmp.reset();
  etmp.append(projPtr->produceIon());
  etmp.append(targPtr->produceIon());
  double mA = projPtr->mN();
  double mB = targPtr->mN();
  double eCM = beamSetupPtr->eCM;
  double pz = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB)
                           * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;

  etmp[1].p(max(projPtr->A(), 1)*Vec4(0.0, 0.0, pz,
                                      sqrt(pow2(pz) + pow2(mA))));
  etmp[1].m(particleDataPtr->m0(idProj));
  etmp[2].p(max(targPtr->A(), 1)*Vec4(0.0, 0.0, -pz,
                                      sqrt(pow2(pz) + pow2(mB))));
  etmp[2].m(particleDataPtr->m0(idTarg));
  etmp[0].p(etmp[1].p() + etmp[2].p());
  etmp[0].m(etmp[0].mCalc());

}

//--------------------------------------------------------------------------

// Take all sub-events and merge them together.

bool Angantyr::buildEvent(list<EventInfo> & subEventsIn) {

  resetEvent();
  Event & etmp = pythia[HADRON]->event;
  double bx = 0.5*FM2MM*hiInfo.b()*cos(hiInfo.phi());
  double by = 0.5*FM2MM*hiInfo.b()*sin(hiInfo.phi());
  etmp[1].vProd( bx,  by, 0.0, 0.0);
  etmp[2].vProd(-bx, -by, 0.0, 0.0);

  // Start with the signal event(s)
  if ( hasSignal ) {
    bool found = false;
    for ( list<EventInfo>::iterator sit = subEventsIn.begin();
          sit != subEventsIn.end(); ++sit  ) {
      // Make sure the events are not SoftQCD or LowEnergyQCD events.
      if ( sit->code >= 101 && sit->code <= 109 ) continue;
      if ( sit->code >= 151 && sit->code <= 159 ) continue;
      addSubEvent(etmp, sit->event);
      hiInfo.select(sit->info);
      subEventsIn.erase(sit);
      found = true;
      break;
    }
    if ( !found ) {
      loggerPtr->ERROR_MSG("failed to generate signal event");
      return false;
    }
  } else
    hiInfo.select(subEventsIn.begin()->info);

  // Then all the others
  for ( list<EventInfo>::iterator sit = subEventsIn.begin();
        sit != subEventsIn.end(); ++sit  ) {
    addSubEvent(etmp, sit->event);
  }
  // Add statistics about participating nucleons and subcollisions.
  hiInfo.glauberStatistics();

  // Finally add all nucleon remnants.
  return addNucleusRemnants();
}

//--------------------------------------------------------------------------

// Construct nucleus remnants from all non-interacting nucleons and
// add them to the main event.

bool Angantyr::addNucleusRemnants() {

  Event & etmp = pythia[HADRON]->event;
  BeamSetup & bs = *beamSetupPtr;
  ParticleData & pdt = pythia[HADRON]->particleData;

  // Get beam particle energies in rest frame.
  double eA = 0.5*(pow2(bs.eCM) + pow2(bs.mA) - pow2(bs.mB))/bs.eCM;
  double eB = bs.eCM - eA;


  // (1) Sum up number of non-interacted nucleons in the
  // projectile and target. Add them directly if they are not protons or
  // neutrons.
  int npp = 0;
  int nnp = 0;
  for (const Nucleon& nucleon : proj)
    if (!nucleon.event()) {
      if ( abs(nucleon.id()) == 2212 ) ++npp;
      else if ( abs(nucleon.id()) == 2112 ) ++nnp;
      else etmp.append(nucleon.id(), 14, 1, 0, 0, 0, 0, 0,
                       0.0, 0.0, sqrt(pow2(eA) - pow2(bs.mA)), eA, bs.mA);
    }

  int npt = 0;
  int nnt = 0;
  for (const Nucleon& nucleon : targ)
    if (!nucleon.event()) {
      if ( abs(nucleon.id()) == 2212 ) ++npt;
      else if ( abs(nucleon.id()) == 2112 ) ++nnt;
      else etmp.append(nucleon.id(), 14, 2, 0, 0, 0, 0, 0,
                       0.0, 0.0, -sqrt(pow2(eB) - pow2(bs.mB)), eB, bs.mB);
    }

  // If every single target and projectile nucleon has interacted,
  // there are no nucleus remanants to add.
  if ( npp + nnp + npt + nnt == 0 ) return true;

  // (2) Sum up the missing momentum: ptot = p_total - sum(p_final).
  Vec4 ptot = etmp[0].p();
  int iPosMax = 0, iNegMax = 0;
  double pPosMax = 0.0, pNegMax = 0.0;
  vector<int> iRemP, iRemT;
  for ( int i = 0, N = etmp.size(); i < N; ++i ) {
    if ( etmp[i].status() <= 0 ) continue;
    ptot -= etmp[i].p();
    // Collect candidate remnant entries.
    if ( etmp[i].status() == 63 ||
         etmp[i].status() == 66 || etmp[i].status() == 14 ) {
      switch ( getBeam(etmp, i) ) {
      case 1: iRemP.push_back(i); break;
      case 2: iRemT.push_back(i); break;
      default: break;
      }
    }
    // Collect fallback particles with the largest p+/p- to be used in
    // the basence of remnants.
    if ( etmp[i].pPos() > pPosMax ) {
      iPosMax = i;
      pPosMax = etmp[i].pPos();
    }
    if ( etmp[i].pNeg() > pNegMax ) {
      iNegMax = i;
      pNegMax = etmp[i].pNeg();
    }
  }
  if ( iRemP.empty() ) iRemP.push_back(iPosMax);
  if ( iRemT.empty() ) iRemT.push_back(iNegMax);

  // (3) Decide whether projectile/target spectators form a nucleus
  // remnants:
  Vec4 pAr, pBr;
  int idAr = 0, idBr = 0;
  double mAr = 0.0, mBr = 0.0;

  if ( npp + nnp > 1 ) {
    // More than one nucleon: we need a composite remnant.
    idAr = 1000000009 + 10000*npp + 10*(nnp + npp);
    //    mAr = (npp + nnp)*bs.mA;
    pdt.addParticle(idAr, "NucRem", 0, 3*npp, 0, (npp + nnp)*bs.mA);
    pdt.particleDataEntryPtr(idAr)->setHasChanged(false);
  }
  // Only one spectator: Simply use p oe n.
  else if ( npp == 1 ) idAr = 2212;
  else if ( nnp == 1 ) idAr = 2112;
  // If no nucleus remnants, use nucleon remnants for energy conservation.
  else {
    for ( int i : iRemP ) pAr += Vec4(0.0, 0.0, etmp[i].pz(), etmp[i].e());
    ptot += pAr;
  }
  if ( bs.idA < 0 ) idAr = -idAr;

  if ( npt + nnt > 1 ) {
    // More than one necleon: we need a composite remnant.
    idBr = 1000000009 + 10000*npt + 10*(nnt + npt);
    pdt.addParticle(idBr, "NucRem", 0, 3*npt, 0, (npt + nnt)*bs.mB);
    pdt.particleDataEntryPtr(idBr)->setHasChanged(false);
  }
  // Only one spectator: Simply use p oe n.
  else if ( npt == 1 ) idBr = 2212;
  else if ( nnt == 1 ) idBr = 2112;
  // If no nucleus remnants, use nucleon remnants energy conservation.
  else {
    for ( int i : iRemT ) pBr += Vec4(0.0, 0.0, etmp[i].pz(), etmp[i].e());
    ptot += pBr;
  }
  if ( bs.idB < 0 ) idBr = -idBr;

  // (4) Do the two body kinematics in the restframe of remnants to be added.
  double mABr = ptot.mCalc();
  auto M = fromCMframe(ptot);
  mAr = idAr ? pdt.m0(idAr) : pAr.mCalc();
  mBr = idBr ? pdt.m0(idBr) : pBr.mCalc();
  double eAr = 0.5*(pow2(mABr) + pow2(mAr) - pow2(mBr))/mABr;
  double eBr = mABr - eAr;
  // Fail if two-body knematics is not possible.
  if ( eAr < mAr || eBr < mBr ) return false;

  // (5) Add new remnants (status 14) or boost existing ones.
  // For the projectile,
  Vec4 pArp = M*Vec4(0.0, 0.0,  sqrt(pow2(eAr) - pow2(mAr)), eAr);
  if ( idAr ) etmp.append(idAr, 14, 1, 0, 0, 0, 0, 0, pArp, mAr);
  else {
    auto MA = fromCMframe(pArp)*toCMframe(pAr);
    for ( int i : iRemP ) etmp[i].rotbst(MA);
  }

  // and for the target.
  Vec4 pBrp = M*Vec4(0.0, 0.0, -sqrt(pow2(eBr) - pow2(mBr)), eBr);
  if ( idBr ) etmp.append(idBr, 14, 2, 0, 0, 0, 0, 0, pBrp, mBr);
  else {
    auto MB = fromCMframe(pBrp)*toCMframe(pBr);
    for ( int i : iRemT ) etmp[i].rotbst(MB);
  }

  // Temporary check for energy-momentum inbalance.
  if ( mode("Angantyr:testMode") ) {
    ptot = etmp[0].p();
    for ( auto p : etmp ) if ( p.isFinal() ) ptot -= p.p();
    if ( ptot.pAbs() + abs(ptot.e()) > 1.0e-6*etmp[0].e() )
      cout << "Energy imbalance error in nucleus remnants " << 1000.0*ptot
           << " recoilers " << iPosMax << "," << iNegMax << endl;
  }

  // Return successful.
  return true;

}

//--------------------------------------------------------------------------

// Set beam kinematics.

bool Angantyr::setKinematics(){
  unifyFrames();
  doLowEnergyNow = false;
  if ( beamSetupPtr->eCM < eCMlow ) {
    if ( !doLowEnergy ) return false;
    doLowEnergyNow = true;
  }
  if ( doLowEnergyNow ) {
    lowEnergyCollPtr->
      updateSig(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM);
    hiInfo.avNDbSave = lowEnergyCollPtr->avNDB();
  } else {
    if ( !sigTotNN.calc(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM))
      return false;
    collPtr->updateSig(projPtr->idN(), targPtr->idN(), beamSetupPtr->eCM);
    if ( !collPtr->setKinematics(beamSetupPtr->eCM) ) return false;
    hiInfo.avNDbSave = collPtr->avNDB();
  }

  // If cascade mode update cross sections for all nuclei in the
  // medium.
  if ( cascadeMode ) updateMedium();

  bGenPtr->updateWidth();
  projPtr->setPN(beamSetupPtr->pAinit);
  targPtr->setPN(beamSetupPtr->pBinit);
  return true;
}

bool Angantyr::setKinematicsCM() {
  hiInfo.glauberReset();
  if ( !setKinematics() ) return false;
  if (!glauberOnly && sabsMode >= 0 &&
      !pythia[SASD]->setKinematics(beamSetupPtr->eCM) )
    return false;
  return pythia[MBIAS]->setKinematics(beamSetupPtr->eCM);
}


bool Angantyr::setKinematics(double eCMIn) {
  if ( eCMIn == beamSetupPtr->eCM ) return true;
  if ( !beamSetupPtr->setKinematics(eCMIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(double eAIn, double eBIn) {
  if ( eAIn == beamSetupPtr->eA && eBIn == beamSetupPtr->eB )
    return true;
  if ( !beamSetupPtr->setKinematics(eAIn, eBIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(double pxAIn, double pyAIn, double pzAIn,
  double pxBIn, double pyBIn, double pzBIn) {
  if ( pxAIn == beamSetupPtr->pxA && pyAIn == beamSetupPtr->pyA &&
       pzAIn == beamSetupPtr->pzA && pxBIn == beamSetupPtr->pxB &&
       pyBIn == beamSetupPtr->pyB && pzBIn == beamSetupPtr->pzB )
    return true;

  if ( !beamSetupPtr->setKinematics(pxAIn, pyAIn, pzAIn,
                                    pxBIn, pyBIn, pzBIn) ) return false;
  return setKinematicsCM();
}

bool Angantyr::setKinematics(Vec4 pA, Vec4 pB) {
  return setKinematics(pA.px(), pA.py(), pA.pz(), pB.px(), pB.py(), pB.pz());
}

//--------------------------------------------------------------------------

// Make sure the correct information is available irrespective of frame type.

void Angantyr::unifyFrames() {
  BeamSetup &bs = *beamSetupPtr;

  if ( bs.frameType == 1 ) {
    bs.eA     = bs.eB = bs.eCM/2;
    bs.pzA    =  sqrt(pow2(bs.eA) - pow2(bs.mA));
    bs.pzB    = -sqrt(pow2(bs.eB) - pow2(bs.mB));
    bs.pxA    = bs.pyA = bs.pxB = bs.pyB = 0.0;
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
  } else if ( bs.frameType == 3 ) {
    bs.eA     = sqrt(pow2(bs.pxA) + pow2(bs.pyA) + pow2(bs.pzA) + pow2(bs.mA));
    bs.eB     = sqrt(pow2(bs.pxB) + pow2(bs.pyB) + pow2(bs.pzB) + pow2(bs.mB));
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
    bs.eCM    = (bs.pAinit + bs.pBinit).mCalc();
  } else {
    // If beam energy is set to less than the mass, it is assumed at rest.
    if (bs.eA < bs.mA ||
        ( projPtr && projPtr->A() > 1 &&
          bs.eA <= particleDataPtr->m0(2112) ) ) {
      bs.pzA = 0.;
      bs.eA = bs.mA;
    }
    else {
      bs.pzA = sqrt(pow2(bs.eA) - pow2(bs.mA));
    }
    if ( bs.eB <= bs.mB ||
         ( targPtr && targPtr->A() > 1 &&
           bs.eB <= particleDataPtr->m0(2112) ) ) {
      bs.pzB = 0.;
      bs.eB = bs.mB;
    }
    else {
      bs.pzB = -sqrt(pow2(bs.eB) - pow2(bs.mB));
    }

    bs.pxA    = bs.pyA = bs.pxB = bs.pyB = 0.0;
    bs.pAinit = Vec4(bs.pxA, bs.pyA, bs.pzA, bs.eA);
    bs.pBinit = Vec4(bs.pxB, bs.pyB, bs.pzB, bs.eB);
    bs.eCM    = (bs.pAinit + bs.pBinit).mCalc();
  }

  if ( !bs.doMomentumSpread ) {
    bs.pAnow = bs.pAinit;
    bs.pBnow = bs.pBinit;
  }

}

//--------------------------------------------------------------------------

// The main method called from Pythia::next().

bool Angantyr::next() {

  if (doSDTest)
    return nextSASD(104);

  pythia[HADRON]->event.clear();
  int itry = MAXTRY;
  bool first = true;

  while ( itry-- && !doAbort) {

    // If cascade mode we only get one try.
    if ( cascadeMode && !first ) return false;
    first = false;

    // Generate impact parameter, nuclei, and sub-collisions.
    double bweight = 0.0;
    Vec4 bvec = bGenPtr->generate(bweight);

    proj = Nucleus(projPtr->generate(), bvec / 2.);
    targ = Nucleus(targPtr->generate(), -bvec / 2.);
    if ( doLowEnergyNow ) {
      subColls = lowEnergyCollPtr->getCollisions(proj, targ);
    } else {
      collPtr->generateNucleonStates(proj, targ);
      subColls = collPtr->getCollisions(proj, targ);
    }

    hiInfo.setSubCollisions(&subColls);
    hiInfo.addAttempt(subColls.T(), bvec.pT(), bvec.phi(),
                      bweight, bGenPtr->xSecScale());

    if ( subColls.empty() ) continue;
    if ( glauberOnly ) return true;

    list<EventInfo> subEvents;

    // Veto early if the primary sub-collision is unwanted.
    auto isVetoed = [&]() -> bool {
      return ( !subEvents.empty() &&
           vetoPrimaryProcess.find(subEvents.begin()->info.code()) !=
               vetoPrimaryProcess.end() );
                    };

    if ( !genAbs(subColls, subEvents) ) {
      loggerPtr->WARNING_MSG("could not setup signal or ND collisions");
      continue;
    }
    if ( hasSignal ) {
      if ( subEvents.empty() ) continue;
    } else {
      if ( isVetoed() ) continue;
    }

    // Collect absorptively wounded nucleons in secondary sub-collisions.
    addSASD(subColls);

    // Collect full low-energy annihilation collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::LANN, 158, Nucleon::ABS, Nucleon::ABS)}) ) {
      loggerPtr->ERROR_MSG("could not setup LANN sub-collision");
      continue;
    }

    // Collect full double diffraction collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::DDE, 105, Nucleon::DIFF, Nucleon::DIFF)}) ) {
      loggerPtr->ERROR_MSG("could not setup DD sub-collision");
      continue;
    }
    if ( isVetoed() ) continue;

    // Collect full low-energy excitation collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::LEXC, 157, Nucleon::DIFF, Nucleon::DIFF)}) ) {
      loggerPtr->ERROR_MSG("could not setup LEXC sub-collision");
      continue;
    }
    if ( isVetoed() ) continue;

    // Collect full single diffraction collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::SDEP, 103, Nucleon::DIFF, Nucleon::ELASTIC),
         CollDesc(SubCollision::SDET, 104, Nucleon::ELASTIC, Nucleon::DIFF)
        }) ) {
      loggerPtr->ERROR_MSG("could not setup SD sub-collision");
      continue;
    }
    if ( isVetoed() ) continue;

    // Collect secondary single diffractive sub-collisions.
    addSecondaries(subColls,
      {CollDesc(SubCollision::DDE,  103, Nucleon::DIFF, Nucleon::UNWOUNDED ),
       CollDesc(SubCollision::SDEP, 103, Nucleon::DIFF, Nucleon::UNWOUNDED ),
       CollDesc(SubCollision::DDE,  104, Nucleon::UNWOUNDED, Nucleon::DIFF ),
       CollDesc(SubCollision::SDET, 104, Nucleon::UNWOUNDED, Nucleon::DIFF )});

    // Collect secondary single diffractive excitations as elastic if
    // the excited one is already excited.
    addSecondaries(subColls,
      {CollDesc(SubCollision::SDEP, 102, Nucleon::UNWOUNDED, Nucleon::ELASTIC),
       CollDesc(SubCollision::SDET, 102, Nucleon::ELASTIC, Nucleon::UNWOUNDED)
      });

    // Collect full central diffraction collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::CDE, 106, Nucleon::ELASTIC, Nucleon::ELASTIC)
        }) ) {
      loggerPtr->ERROR_MSG("could not setup CD sub-collision");
      continue;
    }
    if ( isVetoed() ) continue;

    // Collect full low-energy resonance collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::LRES, 159, Nucleon::ABS, Nucleon::ABS)}) ) {
      loggerPtr->ERROR_MSG("could not setup LRES sub-collision");
      continue;
    }
    if ( isVetoed() ) continue;

    // Collect secondary central diffractive sub-collisions.
    // addCDsecond(subColls);
    addSecondaries(subColls,
      {CollDesc(SubCollision::CDE, 106, Nucleon::ELASTIC, Nucleon::ELASTIC)});

    // Collect full elastic collisions.
    if ( !addSubCollisions(subColls, subEvents,
        {CollDesc(SubCollision::ELASTIC, 102, Nucleon::ELASTIC,
            Nucleon::ELASTIC)}) ) {
      loggerPtr->ERROR_MSG("could not setup elastic sub-collision");
      continue;
    }

    // Collect secondary elastic sub-collisions.
    addSecondaries(subColls,
      {CollDesc(SubCollision::ELASTIC, 102, Nucleon::ELASTIC, Nucleon::ELASTIC)
      });
    if ( isVetoed() ) continue;

    // After all sub-events have been collected, bunch them together.
    if ( subEvents.empty() || !buildEvent(subEvents) ) {
      loggerPtr->ERROR_MSG("failed to build full event");
      continue;
    }

    // Hadronise everything, if requested.
    if (doHadronLevel) {
     if ( HIHooksPtr && HIHooksPtr->canForceHadronLevel() ) {
        if ( !HIHooksPtr->forceHadronLevel(*pythia[HADRON]) ) continue;
      } else {
        if ( !pythia[HADRON]->forceHadronLevel(false) ) continue;
      }
    }

    // Finally, boost to the requested frame and optionally do vertex
    // spreading.
    pythia[HADRON]->event.rotbst(
      fromCMframe(beamSetupPtr->pAnow, beamSetupPtr->pBnow));

    if ( settingsPtr->flag("Beams:allowVertexSpread") ) {
      pythia[HADRON]->getBeamShapePtr()->pick();
      Vec4 vertex = pythia[HADRON]->getBeamShapePtr()->vertex();
      for ( Particle & p : pythia[HADRON]->event ) p.vProdAdd( vertex);
    }

    hiInfo.accept();
    updateInfo();
    return true;
  }

  if (doAbort)
    loggerPtr->ABORT_MSG("Angantyr was aborted due to a critical error");
  else
    loggerPtr->ABORT_MSG("too many attempts to generate a working impact "
      "parameter point", "consider reducing "
                         "HeavyIon:bWidth or HeavyIon:bWidthCut ");
  hiInfo.reject();
  return false;

}

//--------------------------------------------------------------------------

// Parse a given sub-event and streamline it to be easier to handle.

bool Angantyr::streamline(EventInfo& ei) {
  Event & ev = ei.event;
  int & code = ei.code;
  vector<int> & premn = ei.projRems[ei.coll->proj];
  vector<int> & tremn = ei.targRems[ei.coll->targ];
  int & ipqel = ei.projEl[ei.coll->proj];
  int & itqel = ei.targEl[ei.coll->targ];

  // Function for sorting in pseudo rapidity.
  auto etasort =  [](const Particle & a, const Particle & b)
                        { return a.eta() > b.eta(); };

  // Append a particle but modify its status and inheritance.
  auto append =
    [] (Event & e, const Particle & p, int s, int m1, int m2, int d1, int d2) {
      e.append(p);
      e.back().status(s);
      e.back().mothers(m1, m2);
      e.back().daughters(d1, d2);
    };

  // Convert an id into its diffractive counterpart.
  // (eg. 2212 -> 9902210)
  auto difId = [](int id) { int idnew = 10*(abs(id)/10) + 9900000;
    return id > 0? idnew: -idnew;};

  // If we have low energy versions of minimum bias events, remove the
  // hadrons produced in fragmentation, and massage mother/daughter
  // links so that it looks more like the corresponing HardQCD event
  // if possible.
  if ( code == 151 ) {

    // Special case of too low energy.
    if ( ev.size() < 7 ) return true;

    vector<Particle> partons = { ev[3], ev[4], ev[5], ev[6] };
    for ( auto & p : partons ) {
      if ( p.daughter1() == p.daughter2() || !p.colType() ) {
        // Special case of too low energy strings, leave as it is and
        // hope for the best.
        return true;
      }
    }
    ev.popBack(ev.size() - 3);
    ev[1].daughters(3);
    ev[2].daughters(4);
    // Make the most forward incoming parton and outgoing remnant come
    // from the first beam.
    sort(partons.begin(), partons.end(), etasort);
    append(ev, partons[1], -21, 1, 0, 5, 6);
    append(ev, partons[2], -21, 2, 0, 5, 6);
    append(ev, partons[1],  23, 3, 4, 0, 0);
    append(ev, partons[2],  23, 3, 4, 0, 0);
    append(ev, partons[0],  63, 1, 0, 0, 0);
    append(ev, partons[3],  63, 2, 0, 0, 0);
    ipqel = itqel = 0;
    premn = {7};
    tremn = {8};
  } else if ( code == 153 ) {
    Particle pscatp = ev[3];
    Particle premnp = ev[4];
    if ( pscatp.eta() > premnp.eta() ) swap(pscatp, premnp);
    Particle elscat = ev[5];
    ev.popBack(ev.size() - 3);
    ev[1].daughters(3);
    ev[2].daughters(4);
    Vec4 pDif = pscatp.p() + premnp.p();
    ev.append(difId(ev[1].id()), -15, 1, 0, 5, 6, 0, 0, pDif, pDif.mCalc());
    append(ev, elscat, 14, 2, 0, 0, 0);
    append(ev, pscatp,  24, 3, 0, 0, 0);
    append(ev, premnp,  63, 3, 0, 0, 0);
    ipqel = 0;
    itqel = 4;
    premn = {6};
    tremn = {4};
  } else if ( code == 154 ) {
    Particle pscatt = ev[4];
    Particle premnt = ev[5];
    if ( pscatt.eta() < premnt.eta() ) swap(pscatt, premnt);
    Particle elscat = ev[3];
    ev.popBack(ev.size() - 3);
    ev[1].daughters(3);
    ev[2].daughters(4);
    append(ev, elscat, 14, 1, 0, 0, 0);
    Vec4 pDif = pscatt.p() + premnt.p();
    ev.append(difId(ev[2].id()), -15, 2, 0, 5, 6, 0, 0, pDif, pDif.mCalc());
    append(ev, pscatt, 24, 4, 0, 0, 0);
    append(ev, premnt, 63, 4, 0, 0, 0);
    ipqel = 3;
    itqel = 0;
    premn = {3};
    tremn = {6};
  } else if ( code == 155 ) {
    Particle pscatp = ev[3];
    Particle premnp = ev[4];
    if ( pscatp.eta() > premnp.eta() ) swap(pscatp, premnp);
    Particle pscatt = ev[5];
    Particle premnt = ev[6];
    if ( pscatt.eta() < premnt.eta() ) swap(pscatt, premnt);
    ev.popBack(ev.size() - 3);
    ev[1].daughters(3);
    ev[2].daughters(4);
    Vec4 pDifp = pscatp.p() + premnp.p();
    ev.append(difId(ev[1].id()), -15, 1, 0, 5, 6, 0, 0, pDifp, pDifp.mCalc());
    Vec4 pDift = pscatt.p() + premnt.p();
    ev.append(difId(ev[2].id()), -15, 2, 0, 7, 8, 0, 0, pDift, pDift.mCalc());
    append(ev, pscatp, 24, 3, 0, 0, 0);
    append(ev, premnp, 63, 3, 0, 0, 0);
    append(ev, pscatt, 24, 4, 0, 0, 0);
    append(ev, premnt, 63, 4, 0, 0, 0);
    ipqel = itqel = 0;
    premn = {6};
    tremn = {8};
  } else if ( code == 157 || code == 102 || code == 152 ) {
    ev[1].daughters(3);
    ev[2].daughters(4);
    ev[3].mothers(1);
    ev[4].mothers(2);
    ipqel = 3;
    itqel = 4;
    premn = {3};
    tremn = {4};
  } else if ( code == 158 ) {
    ev.popBack(ev.size() - 7);
    while ( ev.back().colType() == 0 ) ev.popBack();
    ipqel = itqel = 0;
    for ( int i = 3; i < ev.size(); ++i ) {
      ev[i].statusPos();
      ev[i].daughters();
      if ( ev[1].id()*ev[i].id() > 0 )
        premn.push_back(i);
      else
        tremn.push_back(i);
    }
  } else if ( code == 159 ) {
    ipqel = itqel = 0;
    premn.clear();
    tremn.clear();
  } else {
    // Special for central diffractive excitation.
    if ( code == 106 ) {
      ev[3].mothers(1);
      ev[4].mothers(2);
    }

    // Special for secondary absorptives.
    if ( code == 101 && ev[3].status() == 14 ) code = 104;
    else if ( code == 101 && ev[4].status() == 14 ) code = 103;

    // Look through the event to find remnants to act as
    // recoilers, if any. Also look for quasi-elastically scattered beam
    // particles.
    for ( int i = 3, N = ev.size(); i < N; ++i ) {
      if ( !ev[i].isFinal() ) continue;
      if ( int irem = isRemnant(ei, i) ) {
        if ( code == 103 && ( irem == 3 || irem == 5 ) ) irem = 1;
        if ( code == 104 && ( irem == 4 || irem == 6 ) ) irem = 2;
        if ( code == 105 && irem > 2 ) {
          if ( ev[irem].status() == -15 ) irem = ev[irem].mother1();
          if ( ev[irem].status() == -13 &&  ev[irem].id() != 990 &&
               ev[ev[irem].mother1()].status() == -15 )
            irem = ev[ev[irem].mother1()].mother1();
        }
        if ( irem == 1 ) premn.push_back(i);
        if ( irem == 2 ) tremn.push_back(i);
      }
      if ( ev[i].mother1() == 1 && ev[i].id() == ev[1].id() ) ipqel = i;
      if ( ev[i].mother1() == 2 && ev[i].id() == ev[2].id() ) itqel = i;
    }
    // (quasi) elastically scattered nucleons are always also remnants.
    if ( ipqel ) premn = {ipqel};
    if ( itqel ) tremn = {itqel};
  }

  return true;
}

//--------------------------------------------------------------------------

// Insert n entries before pos in the event.

void Angantyr::insertEntries(Event & e, int pos, int n) {

  if ( pos >= e.size() || pos < 0 ) return;

  for ( int i = 0; i < n; ++i ) e.append(Particle());
  for ( int i = e.size() - 1; i >= pos + n; --i )
    (e[i] = e[i - n]).offsetHistory(pos - 1, n, pos - 1, n);
  for ( int i = 0; i< pos; ++i )
    e[i].offsetHistory(pos - 1, n, pos - 1, n);

}

//--------------------------------------------------------------------------

// Helper function  to get net flavour content of a hadron or (di)quark.

valarray<int> getNetFlavour(const Particle & p) {
  valarray<int> ret(0, 5);
  int Q = (p.id() > 0? 1: -1);
  int ida = p.idAbs();
  if ( ida == 0  ) return ret;
  if ( ida < 6 ) {
    ret[ida - 1] += Q;
  } else if ( p.isDiquark() ) {
    ret[(ida/1000)%10 - 1] += Q;
    ret[(ida/100)%10 - 1] += Q;
  } else if ( p.particleDataEntry().isMeson() ) {
    int qh = (ida/100)%10;
    if ( qh%2 ) Q = -Q;
    ret[qh - 1] += Q;
    ret[(ida/10)%10 - 1] -= Q;
  } else if ( p.particleDataEntry().isBaryon() ) {
    ret[(ida/1000)%10 - 1] += Q;
    ret[(ida/100)%10 - 1] += Q;
    ret[(ida/10)%10 - 1] += Q;
  }
  return ret;
}

//--------------------------------------------------------------------------

// Fix up a diffractive event with only one gluon extracted from one
// side to look lie an single diffractive event.

bool Angantyr::fixSecondaryAbsorptive(Event & ev, double xpom) {

  // First decide by the sign of xpom which beam should be elastically
  // scattered, and which should be diffractively excited (obeam).
  bool targetside = ( xpom < 0.0 );
  int beam = targetside? 2: 1;
  int obeam = 3 - beam;

  // We need to keep track of all flavours in the beam to be elastic.
  valarray<int> beamquarks = getNetFlavour(ev[beam]);

  // First trace all colours in the final state and find the two
  // remnants (quark and diquark (or antiquark in the case of a
  // meson)). Also sum up all other momenta.

  // map from colour lines to the indices of the final state endpoints
  // (note special code for junctions below).
  struct ColEndpoints {
    int colIdx  = 0; // particle index or encoded junction ref (<0)
    int acolIdx = 0;
  };
  map<int, ColEndpoints> cols;

  // List of all coloured and anticoloured remnants.
  vector<int> colrems;
  vector<int> acolrems;
  // List of all remnant colours and anticolours of remnants.
  vector<int> remscol;
  vector<int> remsacol;

  // Map colour lines to junction legs.
  for ( int j = 0; j < ev.sizeJunction(); ++j )
    for ( int l = 0; l < 3; ++l )
      if ( ev.kindJunction(j)%2 )
        cols[ev.colJunction(j, l)].acolIdx = -1000*(j + 1) - l;
      else
        cols[ev.colJunction(j, l)].colIdx = -1000*(j +1) - l;

  // Here we collect the momentum of everything that is not the
  // remnants to become an elastically scattered beam.
  Vec4 pRec;

  for ( int ip = 1, Np = ev.size(); ip < Np; ++ip ) {
    auto & p = ev[ip];
    if ( !p.isFinal() ) continue;
    if ( p.status() == 63 && p.mother1() == beam ) {

      // Make sure we find all beam flavours.
      beamquarks -= getNetFlavour(p);

      // Save the colours and indices of the interesting remnants.
      if ( p.col() ) {
        colrems.push_back(ip);
        remscol.push_back(p.col());
      }
      if ( p.acol() ) {
        acolrems.push_back(ip);
        remsacol.push_back(p.acol());
      }
    } else {
      // Sum momentum of everyting else.
      pRec += p.p();
    }

    // Collect information of all colour lines.
    if ( p.col() ) cols[p.col()].colIdx = ip;
    if ( p.acol() ) cols[p.acol()].acolIdx = ip;
  }

  // Make sure we have found all flavours in the beam, and only
  // remnats consistent with only having extracted gluons from the
  // beam to be elastically scattered.
  if ( beamquarks.min() != 0 || beamquarks.max() != 0 ||
       colrems.size() != 1 || acolrems.size() != 1 ) {
    //    for ( int i = 0; i < 5; ++i ) cout << beamquarks[i];
    loggerPtr->WARNING_MSG("Could not add secondary absorptive");
    return false;
  }

  // We now only have one collered and one anti-coloured remnant.
  int iq   = colrems[0];
  int iqq  = acolrems[0];
  int iqc  = remscol[0];
  int iqqc = remsacol[0];

  int ig = 0;
  // Now go through all cases and make the remnants colour
  // singlet. Optionally keep a gluon, taking some of the momentum.
  if ( iqc == iqqc ) {
    // Remnants in colour singlet - nothing to do. Currently this will
    // never happen, but may be useful in the future.
  } else if ( sabsMode < -1 ) {
    // Change the coloured remnat to a gluon to take some recoil.
    ig = iq;
    ev[ig].id(21);
    ev[ig].acol(ev[iqq].acol());
  } else if ( cols[iqqc].colIdx == cols[iqc].acolIdx ) {
    // Remnants connected to a single gluon, we need to connect that
    // elsewhere. Insert in a dipole so that the (Ariadne definition
    // of the) transverse momentum is as small as possible.
    int igi = cols[iqqc].colIdx;
    const Vec4 & p2 = ev[igi].p();
    double ptmin = ev[0].m();
    int colsel = 0;
    for ( auto it : cols ) {
      if ( it.first == iqqc || it.first == iqc ) continue;
      if ( cols[it.first].colIdx < 0 || cols[it.first].acolIdx < 0 ) continue;
      const Vec4 & p1 = ev[it.second.colIdx].p();
      const Vec4 & p3 = ev[it.second.acolIdx].p();
      double pt = (p1 + p2).mCalc()*(p2 + p3).mCalc()/(p1 + p2 + p3).mCalc();
      if ( pt < ptmin ) {
        colsel = it.first;
        ptmin = pt;
      }
    }
    if ( colsel == 0 ) return false;
    ev[igi].col(colsel);
    ev[cols[colsel].colIdx].col(ev[igi].acol());
    ev[iq].col(ev[iqq].acol());
  } else if ( cols[iqc].colIdx < 0 && cols[iqqc].acolIdx < 0 ) {
    // Both remnants connected to junctions won't work.
    return false;
  } else if ( cols[iqqc].colIdx < 0 ) {
    // Remnants connecting to junction is tricky.
    Junction & J = ev.getJunction((-cols[iqqc].colIdx)/1000 - 1);
    // Should never happen.
    if ( J.kind()%2 != 0 ) return false;
    int l = (-cols[iqqc].colIdx)%1000;
    J.col(l, iqc);
    ev[iq].col(iqqc);
  } else if ( cols[iqc].acolIdx < 0 ) {
    // Remnants connecting to junction is tricky.
    int j = (-cols[iqc].acolIdx)/1000 - 1;
    Junction & J = ev.getJunction(j);
    // Should never happen.
    if ( J.kind()%2 != 1 ) return false;
    int l = (-cols[iqc].acolIdx)%1000;
    J.col(l, iqqc);
    ev[iqq].acol(iqc);
  } else {
    // Remnants connected normally. We just need to reconnect one.
    ev[cols[iqc].acolIdx].acol(iqqc);
    ev[iqq].acol(iqc);
  }

  // Now replace the remnants with a proton. Fix up the colours and
  // momenta, and set up the event to be a single diffractive one. We
  // calculate:
  // the original summed remnant momentum,
  Vec4 pRem = ev[iq].p() + ev[iqq].p();
  // the momentum of the pomeron to be inserted.
  Vec4 ppom;
  // The momenum of the elastically scattered beam,
  Vec4 pqel;
  // and the boost needed for all other particles to get
  // energy-momentum conservarion.
  RotBstMatrix MRec;

  // Optionally we include a gluon so that the mass of the diffractive
  // system corresponds to the exchange of a pomeron with energy
  // fraction xpom.
  if ( ig ) {
    int beamdir = ( targetside? -1: 1 );
    ppom.e(abs(xpom)*ev[beam].e());
    pqel.e(ev[beam].e() - ppom.e());
    double z2 = pow2(pqel.e()) - ev[beam].m2Calc();
    if ( z2 < 0.0 ) return false;
    pqel.pz(beamdir*sqrt(z2));
    ppom.pz(ev[beam].pz() - pqel.pz());
    Vec4 p1move(0.0, 0.0, ppom.pz() + ev[obeam].pz() - pRec.pz(),
                ppom.e() + ev[obeam].e() - pRec.e());
    Vec4 p2move(0.0, 0.0, pRec.pz(), pRec.e());
    if ( !pShift(p1move, p2move, pRem.pT(), sqrt(pRec.pPos()*pRec.pNeg())) )
      return false;

    ev[ig].p(-pRec.px(), -pRec.py(), p1move.pz(), p1move.e());
    ev[ig].id(21);
    ev[ig].m(0);
    MRec.bst(0.0, 0.0, -pRec.pz()/pRec.e());
    MRec.bst(0.0, 0.0, p2move.pz()/p2move.e());

    ev.remove(iqq, iqq);
    if ( ig > iqq ) --ig;

  } else {

    Vec4 p1move(0.0, 0.0, pRem.pz(), pRem.e());
    Vec4 p2move(0.0, 0.0, pRec.pz(), pRec.e());
    if ( !pShift(p1move, p2move, sqrt(pRem.pT2() + pow2(ev[beam].m())),
                 sqrt(pRec.pT2() + pRec.m2Calc())) )
      return false;
    p1move.px(pRem.px());
    p1move.py(pRem.py());
    pqel = p1move;
    ppom = ev[beam].p() - p1move;
    MRec.bst(0.0, 0.0, -pRec.pz()/pRec.e());
    MRec.bst(0.0, 0.0, p2move.pz()/p2move.e());
    ev.remove(max(iq,iqq), max(iq,iqq));
    ev.remove(min(iq,iqq), min(iq,iqq));
  }

  // Finally we manipulate the Event to make it look like a single
  // diffractive excitation event.
  Vec4 pdif = ev[obeam].p() + ppom;

  insertEntries(ev, 3, 4);
  int iqel = beam + 2;
  int idif = 5 - beam;
  int ipom = beam + 4;
  int iocp = 7 - beam;
  ev[iqel] = Particle(ev[beam].id(), 14, beam, 0, 0, 0, 0, 0,
                      pqel, ev[beam].m());
  ev.setEvtPtr(iqel);

  int iddif = (990000 + ev[obeam].idAbs()/10)*10*(ev[obeam].id() < 0? -1: 1);
  ev[idif] = Particle(iddif, -15, obeam, 0, 0, 0, 0, 0,
                      pdif, pdif.mCalc());
  ev.setEvtPtr(idif);

  ev[ipom] = Particle(990, -13, idif, 0, ev[beam].daughter1(), 0, 0, 0,
                      ppom, ppom.mCalc());
  ev.setEvtPtr(ipom);

  ev[iocp] = Particle(ev[obeam].id(), -13, idif, 0, ev[obeam].daughter1(),
                      0, 0, 0,
                      ev[obeam].p(), ev[obeam].m());
  ev.setEvtPtr(iocp);

  ev[beam].daughters(iqel);
  ev[obeam].daughters(idif);
  ev[idif].daughters(ipom, iocp);

  for ( int ip = 7; ip < ev.size(); ++ip ) {
    if ( ip != ig + 4 ) ev[ip].rotbst(MRec);
    if ( ev[ip].mother1() == beam ) ev[ip].mother1(ipom);
    if ( ev[ip].mother2() == beam ) ev[ip].mother2(ipom);
    if ( ev[ip].mother1() == obeam ) ev[ip].mother1(iocp);
    if ( ev[ip].mother2() == obeam ) ev[ip].mother2(iocp);
  }

  return true;

}

//--------------------------------------------------------------------------

// Function used in ProcessSelectorHook to veto processes.

bool Angantyr::ProcessSelectorHook::procveto(const Event& ev) {

  // If primary scattering check x incoming and process number,
  if ( proc > 0 && infoPtr->code()%10 != proc%10 ) return true;
  giveUp = false;
  hasQuark = false;
  iAllowQuarkTries = 0;
  if ( proc%10 == 1 ) return xveto(2, ev, 3, 0);
  if ( proc%10 == 3 && xmax < 0.0 &&
       ev[4].pNeg()/ev[0].pNeg() < 1.0 + xmax ) return true;
  if ( proc%10 == 4 && xmax > 0.0 &&
       ev[3].pPos()/ev[0].pPos() < 1.0 - xmax ) return true;
  return false;

}

//--------------------------------------------------------------------------

// Function used in ProcessSelectorHook to veto final parton level.

bool Angantyr::ProcessSelectorHook::finalveto(const Event& ev) {
  if (xmax == 0.0 || proc%10 != 1) return false;
  int ib = ( xmax > 0? 1: 2);
  bool looksok = true;
  for ( auto & p : ev ) {
    if ( p.status() < 0 && p.mother1() == ib && p.id() != 21 ) {
      looksok = false;
      break;
    }
  }
  if ( ( looksok && ( giveUp || hasQuark ) ) ||
       ( !looksok && !giveUp && !hasQuark ) ) {
    loggerPtr->ERROR_MSG("Something went wrong when generating restricted ND");
    return true;
  }

  return giveUp || hasQuark;
}

//--------------------------------------------------------------------------

// Function used in ProcessSelectorHook to veto emissions.

bool Angantyr::ProcessSelectorHook::
xveto(int type, const Event & event, int oldsize, int iSys) {

  // No veto if no limit has been specified
  if ( xmax == 0.0 || proc%10 != 1 ) return checkVeto(false);

  // Give up if we cannot use this event.
  if ( giveUp ) return checkVeto(true);

  // Scan (the new part of) the event recors to find chaged
  // incoming partons and their momenta and check if they are
  // gluons.
  Vec4 pIn = {};
  set<int> changed;
  int inew = 0;
  for ( int i = oldsize; i < event.size(); ++i ) {
    if ( event[i].status() == -31 ||
         event[i].status() == -21 ) {
      pIn += event[i].p();
      if ( event[i].pz()*xmax > 0.0 ) inew = i;
    }
    if ( event[i].status() == -41 || event[i].status() == -42 ||
         event[i].status() == -53 || event[i].status() == -54 ) {
      pIn += event[i].p();
      changed.insert(event[i].daughter2());
    }
    if ( event[i].status() == -41 && event[i].pz()*xmax > 0.0 ) inew = i;
  }

  // Scan through the parton systems to find the incoming momenta
  // and flavours.
  if ( type != 2 )
    for ( int is = 0; is < partonSystemsPtr->sizeSys(); ++is ) {
      if ( int ia = partonSystemsPtr->getInA(is) )
        if ( changed.find(ia) == changed.end() ) pIn += event[ia].p();
      if ( int ib = partonSystemsPtr->getInB(is) )
        if ( changed.find(ib) == changed.end() ) pIn += event[ib].p();
    }
  if ( mode("Angantyr:testMode") > 1 ) {
    // *** TODO *** This is just to check if we got the correct
    // incoming momenta in the code above. Should be removed.
    Vec4 pf = {};
    for ( auto & p : event ) if ( p.isFinal() ) pf += p.p();

    if ( abs(pIn.pPos() - pf.pPos())/event[0].pPos() > 1.0e-6 ||
         abs(pIn.pNeg() - pf.pNeg())/event[0].pNeg() > 1.0e-6 ) {
      cout << "Angantyr did not get  the exact x-value ina SASD treatment "
           << (event[event.size() - 1].status())/10 << endl;
      cout << setprecision(10) << pIn.pPos()/event[0].pPos() << ", "
           << pIn.pNeg()/event[0].pNeg() << endl;
      cout << pf.pPos()/event[0].pPos() << ", "
           << pf.pNeg()/event[0].pNeg() << endl << endl;
    }
  }

  // Veto if too much momentum has been extracted from the relevant
  if ( ( xmax > 0? pIn.pPos()/event[0].pPos():
         pIn.pNeg()/event[0].pNeg()) > abs(xmax) ) return checkVeto(true);

  // Primary scattering may extract quark, but it needs to evolve
  // back to a gluon in the end.
  if ( type == 2 && event[inew].id() != 21 ) {
    if ( ++iAllowQuarkTries > nAllowQuarkTries ) return checkVeto(true);
    hasQuark = true;
  }

  // Veto if quarks were extracted on the relevant side in
  // secondary scatterings.
  if ( type == 3  && event[inew].id() != 21 ) return checkVeto(true);

  // ISR may evolve a quark to a quark, but only a maximum number of times.
  if ( type == 4  ) {
    if ( event[inew].id() != 21 ) {
      if ( iSys != 0 || !hasQuark ) return checkVeto(true);
      if ( ++iAllowQuarkTries > nAllowQuarkTries )
        return checkVeto(giveUp = true);
    } else
      if ( iSys == 0 ) hasQuark = checkVeto(false);
  }

  return checkVeto(false);

}

//==========================================================================

} // end namespace Pythia8
