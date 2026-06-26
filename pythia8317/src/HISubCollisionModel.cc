// HISubCollisionModel.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the HISubCollisionModel.h header) for
// the built-in heavy ion subcollision models.

#include "Pythia8/Pythia.h"
#include "Pythia8/HISubCollisionModel.h"

namespace Pythia8 {

//==========================================================================

// ImpactParameterGenerator samples the impact parameter space.

//--------------------------------------------------------------------------

// Initialise base class, passing pointers to important objects.

void ImpactParameterGenerator::initPtr(Info & infoIn,
  SubCollisionModel & collIn, NucleusModel & projIn, NucleusModel & targIn) {
  infoPtr = &infoIn;
  settingsPtr = infoIn.settingsPtr;
  rndmPtr = infoIn.rndmPtr;
  loggerPtr = infoIn.loggerPtr;
  collPtr = &collIn;
  projPtr = &projIn;
  targPtr = &targIn;

}

//--------------------------------------------------------------------------

// Initialise base class, may be overridden by subclasses.

bool ImpactParameterGenerator::init() {
  // The width parameter is given in units of femtometer.
  if ( settingsPtr->isParm("HI:bWidth") )
    widthSave = settingsPtr->parm("HI:bWidth");
  else
    widthSave = settingsPtr->parm("HeavyIon:bWidth");

  if ( widthSave <= 0.0 )
    updateWidth();

  cut = settingsPtr->parm("HeavyIon:bWidthCut");
  forceUnitWeight = settingsPtr->flag("HeavyIon:forceUnitWeight");

  // For backward compatibility
  if ( cut == settingsPtr->parmDefault("HeavyIon:bWidthCut") )
    cut = sqrt(-2.0*log(0.01));

  return true;
}

//--------------------------------------------------------------------------

// Set width based on the associated subcollision and nucleus models.

void ImpactParameterGenerator::updateWidth() {
  double Rp = sqrt(collPtr->sigTot()/M_PI)/2.0;
  double RA = max(Rp, projPtr->R());
  double RB = max(Rp, targPtr->R());
  widthSave = RA + RB + 2.0*Rp;
}

//--------------------------------------------------------------------------

// Generate an impact parameter according to a gaussian distribution.

Vec4 ImpactParameterGenerator::generate(double & weight) const {

  if ( forceUnitWeight ) {
    double b = width()*cut*sqrt(rndmPtr->flat());
    double phi = 2.0*M_PI*rndmPtr->flat();
    weight = 1.0;
    return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
  }

  double R = 0.0;
  double Rcut = exp(-pow2(cut)/2.0);
  do
    R = rndmPtr->flat();
  while (R < Rcut);
  double b = sqrt(-2.0*log(R))*width();
  double phi = 2.0*M_PI*rndmPtr->flat();
  weight = 1.0/R;
  return Vec4(b*sin(phi), b*cos(phi), 0.0, 0.0);
}

//==========================================================================

// The SubCollisionModel base class for modeling the collision between
// two nucleons to tell which type of collision has occurred. The
// model may manipulate the corresponding state of the nucleons.

//--------------------------------------------------------------------------

shared_ptr<SubCollisionModel> SubCollisionModel::create(int model) {
  switch (model) {
    case 0: return make_shared<NaiveSubCollisionModel>();
    case 1: return make_shared<DoubleStrikmanSubCollisionModel>();
    case 2: return make_shared<DoubleStrikmanSubCollisionModel>(1);
    case 3: return make_shared<BlackSubCollisionModel>();
    case 4: return make_shared<LogNormalSubCollisionModel>();
    case 5: return make_shared<LogNormalSubCollisionModel>(1);
    case 6: return make_shared<DoubleStrikmanSubCollisionModel>(2);
    default: return nullptr;
  }
}

//--------------------------------------------------------------------------

// Initialize the base class. Subclasses should consider calling this
// in overriding functions.

bool SubCollisionModel::init(int idAIn, int idBIn, double eCMIn) {

  // Store input.
  idASave = idAIn;
  idBSave = idBIn;
  eSave   = eCMIn;

  // Read basic settings.
  NInt = settingsPtr->mode("HeavyIon:SigFitNInt");
  NPop = settingsPtr->mode("HeavyIon:SigFitNPop");
  sigErr = settingsPtr->pvec("HeavyIon:SigFitErr");
  sigFuzz = settingsPtr->parm("HeavyIon:SigFitFuzz");
  fitPrint = settingsPtr->flag("HeavyIon:SigFitPrint");
  impactFudge = settingsPtr->parm("Angantyr:impactFudge");
  elasticMode  = settingsPtr->mode("Angantyr:elasticMode");
  elasticFudge = settingsPtr->parm("Angantyr:elasticFudge");
  eCMlow = settingsPtr->parm("HeavyIon:eCMLowEnergy");

  doVarECM = settingsPtr->flag("Beams:allowVariableEnergy");
  doVarBeams = settingsPtr->flag("Beams:allowIDASwitch");
  if (doVarBeams) {
    idAList = settingsPtr->mvec("Beams:idAList");
    if (idAList.size() == 0) {
      loggerPtr->ABORT_MSG(
        "requested variable beams, but Beams:idAList is empty");
      return false;
    }
    else if (idAList.size() == 1) {
      loggerPtr->WARNING_MSG("requested variable beams, but "
        "Beams:idAList contains only a single entry");
    }
    bool idAIsGood = false;
    for (int idA : idAList) if (idA == idAIn) {
      idAIsGood = true;
      break;
    }
    if (!idAIsGood) {
      loggerPtr->WARNING_MSG("Beams:idA not found in Beams:idAList",
        "defaulting to " + to_string(idAList[0]));
      idASave = idAList[0];
    }
  } else idAList = vector<int>{ idASave };

  if (doVarECM) {
    eMin = settingsPtr->parm("HeavyIon:varECMMin");
    eMax = settingsPtr->parm("HeavyIon:varECMMax");
    eCMPts = settingsPtr->mode("HeavyIon:varECMSigFitNPts");
    if (eMax == 0)
      eMax = eCMIn;
    else if (eMax < eCMIn) {
      loggerPtr->ERROR_MSG("maximum energy is lower than requested eCM");
      return false;
    }
  }
  else {
    eCMPts = 1;
    eMin = eMax = eCMIn;
  }
  updateSig(idAIn, idBIn, eCMIn);

  // If there are parameters, no further initialization is necessary.
  if (nParms() == 0) return true;

  // First try to load configuration from file, if requested.
  int    reuseInitMode = settingsPtr->mode("HeavyIon:SigFitReuseInit");
  string reuseInitFile = settingsPtr->word("HeavyIon:SigFitInitFile");
  bool   reuseWorked   = (reuseInitMode < 0 ||reuseInitMode == 2 ||
                          reuseInitMode == 3) && loadParms(reuseInitFile);

  if (!reuseWorked) {
    if (reuseInitMode == 2) {
      loggerPtr->ABORT_MSG("unable to load parameter data");
      return false;
    }

    // If parameters were not loaded, generate from scratch.
    if (!genParms()) {
      loggerPtr->ABORT_MSG("evolutionary algorithm failed");
      return false;
    }
  }

  // Set parameters at the correct kinematics.
  setKinematics(eCMIn);

  // Set initial avNDb
  if ( settingsPtr->mode("Angantyr:impactMode") < 0.0 || !reuseWorked ) {
    SigEst se = getSig();
    avNDb = se.avNDb;
    avNDolap = se.avNOF;
  }

  // Save parameters to disk, if requested.
  if (reuseInitMode < 0 ||reuseInitMode == 1 ||
      (reuseInitMode == 3 && !reuseWorked) ) {
    if (saveParms(reuseInitFile)) {
      if ( reuseInitMode != -1 ) loggerPtr->INFO_MSG(
        "wrote initialization configuration to file", reuseInitFile);
    }
    else loggerPtr->WARNING_MSG("couldn't save initialization configuration");
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Generate parameters based on run settings and the evolutionary algorithm.

bool SubCollisionModel::genParms() {

  // Initialize with default parameters.
  int nGen = settingsPtr->mode("HeavyIon:SigFitNGen");
  vector<double> defPar = settingsPtr->pvec("HeavyIon:SigFitDefPar");
  if ( settingsPtr->isPVec("HI:SigFitDefPar") )
    defPar = settingsPtr->pvec("HI:SigFitDefPar");
  if (defPar.size() == 0)
    defPar = defParm();
  if (int(defPar.size()) < nParms()) {
    loggerPtr->ERROR_MSG("too few parameters have been specified",
      "(expected " + to_string(nParms())
      + ", got " + to_string(defPar.size()) + ")");
    return false;
  }
  if (int(defPar.size()) > nParms()) {
    loggerPtr->WARNING_MSG("too many parameters have been specified",
      "(expected " + to_string(nParms())
      + ", got " + to_string(defPar.size()) + ")");
    defPar.resize(nParms());
  }

  // Read settings for varECM evolution.
  int doStepwiseEvolve = settingsPtr->mode("HeavyIon:varECMStepwiseEvolve");


  setParm(defPar);

  double sigPPrefND = 0.0;
  if ( doStepwiseEvolve == 2 ) {
    sigTotPtr->calc(2212, 2212, eMax);
    updateSig(2212, idBSave, eMax);
    sigPPrefND = sigTarg[1];
  }

  for (int idANow : idAList) {

    sigTotPtr->calc(idANow, 2212, eMax);
    updateSig(idANow, idBSave, eMax);

    // If doing step-wise evolution in eCM we still need to reset when
    // changing beam type.
    vector<double> parmsNow = defPar;
    if ( doStepwiseEvolve ) {
      // Optionally scale with the non-diffractie cross section.
      if ( doStepwiseEvolve == 2 )
        parmsNow[0] *= sigTarg[1]/sigPPrefND;
      setParm(parmsNow);
    }
    vector<LogInterpolator> subCollParmsNow;

    // If nGen is zero, there is nothing to do, just use the default
    // parameters.
    if (nGen == 0) {
      subCollParmsNow = vector<LogInterpolator>(nParms() + 1);
      for (int iParm = 0; iParm < nParms(); ++iParm)
        subCollParmsNow[iParm] = LogInterpolator(eMin, eMax, {defPar[iParm]});
      subCollParmsNow[nParms()] = LogInterpolator(eMin, eMax, {avNDb});
      subCollParmsMap[idANow] = subCollParmsNow;
      continue;
    }

    // Run evolutionary algorithm.
    if ( fitPrint ) {
      cout << " *------ HeavyIon fitting of SubCollisionModel to "
          << "cross sections ------* " << endl;
      flush(cout);
    }
    if (!evolve(nGen, eMax, idANow)) {
      loggerPtr->ERROR_MSG("evolutionary algorithm failed");
      return false;
    }
    parmsNow = getParm();

    // If we don't care about varECM, we are done.
    if (!doVarECM) {
      if (fitPrint) {
        cout << " *--- End HeavyIon fitting of parameters in "
          << "nucleon collision model ---* "
          << endl << endl;
        cout << " Angantyr Info: To avoid refitting,"
          " add the following lines to your configuration file: " << endl;
        cout << "                HeavyIon:SigFitNGen = 0" << endl;
        cout << "                HeavyIon:SigFitDefAvNDb = " << avNDb << endl;
        cout << "                HeavyIon:SigFitDefPar = ";
        for (int iParm = 0; iParm < nParms(); ++iParm) {
          if (iParm > 0) cout << ",";
          cout << parmsNow[iParm];
        }
        cout << endl << endl;
      }
      subCollParmsNow = vector<LogInterpolator>(nParms() + 1);
      for (int iParm = 0; iParm < nParms(); ++iParm)
        subCollParmsNow[iParm] = LogInterpolator(
          eMin, eMax, {parmsNow[iParm]});
      subCollParmsNow[nParms()] = LogInterpolator(
        eMin, eMax, {avNDb});
      subCollParmsMap[idANow] = subCollParmsNow;
      continue;
    }

    // Vector of size nParms, each entry contains the parameter values.
    vector<vector<double>> parmsByECM(nParms() + 1, vector<double>(eCMPts));

    // Write parameters at original eCM.
    for (int iParm = 0; iParm < nParms(); ++iParm)
      parmsByECM[iParm].back() = parmsNow[iParm];
    // Also store the average non-diffractive impact parameter.
    parmsByECM[nParms()].back() = avNDb;

    // Evolve down to eMin.
    vector<double> eCMs = logSpace(eCMPts, eMin, eMax);
    double eParmScale = pow(eCMs[0]/eCMs[1], 0.12);
    for (int i = eCMPts - 2; i >= 0; --i) {
      // Update to correct eCM.
      double eNow = eCMs[i];
      sigTotPtr->calc(idANow, idBSave, eNow);
      updateSig(idANow, idBSave, eNow);

      // Alternatively reset to default parameters (mostly for debug purposes).
      if (!doStepwiseEvolve)
        setParm(defPar);
      else if (doStepwiseEvolve == 2) {
        for ( double & p : parmsNow ) p *= eParmScale;
        setParm(parmsNow);
      }

      // Evolve and get next set of parameters.
      if (fitPrint)
        cout << " *------------------------------------------"
              "---------------------------* "
            << endl;

      if (!evolve(nGen, eNow, idANow)) {
        loggerPtr->ERROR_MSG("evolutionary algorithm failed");
        return false;
      }
      parmsNow = getParm();
      for (int iParm = 0; iParm < nParms(); ++iParm)
        parmsByECM[iParm][i] = parmsNow[iParm];
      parmsByECM[nParms()][i] = avNDb;

    }
    if (fitPrint){
      cout << " *--- End HeavyIon fitting of parameters in "
          << "nucleon collision model ---* "
          << endl << endl;
      cout << " Angantyr Info: To avoid refitting, you may use the "
              "HeavyIon:SigFitReuseInit parameter "
              "\n                to store the configuration to disk."
           << endl << endl;
    }
    // Reset cross section and parameters to their eCM values.
    sigTotPtr->calc(idASave, idBSave, eMax);
    updateSig(idASave, idBSave, eMax);
    setParm(parmsNow);

    // Store parameter values as logarithmic interpolators.
    subCollParmsNow = vector<LogInterpolator>(nParms() + 1);
    for (int iParm = 0; iParm < nParms() + 1; ++iParm) {
      subCollParmsNow[iParm] = LogInterpolator(eMin, eMax, parmsByECM[iParm]);
      if ( doStepwiseEvolve && iParm < nParms() )
        parmSave[iParm] = parmsByECM[iParm].back();
    }
    subCollParmsMap[idANow] = subCollParmsNow;
  }

  // Set default parameters.
  subCollParmsPtr = &subCollParmsMap.at(idASave);
  for (int iParm = 0; iParm < nParms(); ++iParm)
    parmSave[iParm] = subCollParmsPtr->at(iParm).data().back();
  avNDb = subCollParmsPtr->at(nParms()).data().back();

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Save parameter configuration to settings/disk.

bool SubCollisionModel::saveParms(string fileName) const {

  if (nParms() == 0) {
    loggerPtr->WARNING_MSG("model does not have any parameters");
    return true;
  }

  vector<string> setting;
  ostringstream os;
  os << eCMPts << " " << eMin << " " << eMax << " " << 1;
  setting.push_back(os.str());

  for (int idANow : idAList) {

    // Write idA.
    setting.push_back(to_string(idANow));

    // Each line corresponds to one parameter.
    auto& subCollParmsNow = subCollParmsMap.at(idANow);
    for (int iParm = 0; iParm < nParms() + 1 ; ++iParm) {
      ostringstream oss;
      oss << setprecision(6);
      for (double val : subCollParmsNow[iParm].data())
        oss << " " << val;
      setting.push_back(trimString(oss.str()));
    }
  }

  settingsPtr->wvec("Init:reuseHeavyIonSigFit", setting);

  if ( fileName.length() == 0 ||
       settingsPtr->mode("HeavyIon:SigFitReuseInit") == -1) return true;

  ofstream ofs(fileName);
  if (!ofs.good()) {
    loggerPtr->ERROR_MSG("unable to open file for writing", fileName);
    return false;
  }

  ofs << "Init:reuseHeavyIonSigFit = { " << setting[0] << "," << endl;
  for ( unsigned int i = 1; i < setting.size(); ++i ) {
    ofs << "      " << setting[i];
    if ( i == setting.size() - 1 ) ofs << " }" << endl;
    else ofs << "," << endl;
  }

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Load parameter configuration from settings/disk.

bool SubCollisionModel::loadParms(string fileName) {

  if (nParms() == 0) {
    loggerPtr->WARNING_MSG("model does not have any parameters");
    return true;
  }

  if ( fileName.length() > 0 ) {
    ifstream istest(fileName);
    if ( istest.good() )
      settingsPtr->readFile(fileName);
  }

  vector<string> lines = settingsPtr->wvec("Init:reuseHeavyIonSigFit");
  if ( lines.size() < 2 ) {
    loggerPtr->WARNING_MSG("stored values do not cover requested energy range."
                           " Regenerating.");
    return false;
  }

  // Lambda function as a shorthand for error message.
  auto formatError = [this]() {
    loggerPtr->ERROR_MSG("invalid format");
    return false;
  };

  // Read first line to get energy range and number of interpolation points.
  istringstream is(lines[0]);
  double eMinNow, eMaxNow;
  if ( !( is >> eCMPts >> eMinNow >> eMaxNow) )
    return formatError();
  if (!(eCMPts >= 1) || eMin < 0.999*eMinNow || eMax > eMaxNow*1.001) {
    loggerPtr->ERROR_MSG("stored file does not cover requested energy range");
    return false;
  }
  int fileversion = 0;
  if ( ! ( is >> fileversion ) ) {
    fileversion = 0;
    loggerPtr->WARNING_MSG("reading file with old deprecated format");
  }

  eMin = eMinNow;
  eMax = eMaxNow;

  for ( unsigned int i = 1; i < lines.size();  ) {
    istringstream iss(lines[i++]);
    // Read idA.
    int idANow;
    if  ( !(iss >> idANow) ) return formatError();

    // Read each line and use the data to define an interpolator.
    vector<LogInterpolator> subCollParmsNow(nParms() + 1);
    for (int iParm = 0; iParm < nParms() + 1; ++iParm) {
      istringstream lineStream(lines[i++]);
      vector<double> parmData(eCMPts);
      for (int iPt = 0; iPt < eCMPts; ++iPt) {
        if (!(lineStream >> parmData[iPt])) return formatError();
        // ***TODO *** this should be fixed!
        if ( !fileversion && iParm == nParms() )
          parmData[iPt] /= settingsPtr->parmDefault("Angantyr:impactFudge");
      }

      subCollParmsNow[iParm] = LogInterpolator(eMin, eMax, parmData);
    }

    subCollParmsMap.emplace(idANow, subCollParmsNow);
  }

  // Validate that requested ids have been loaded.
  for (int idANow : idAList) {
    if (subCollParmsMap.find(idANow) == subCollParmsMap.end()) {
      loggerPtr->ERROR_MSG("requested ids not found in stored file");
      return false;
    }
  }

  // Set default parameters.
  subCollParmsPtr = &subCollParmsMap[idASave];
  for (int iParm = 0; iParm < nParms(); ++iParm)
    parmSave[iParm] = subCollParmsPtr->at(iParm).data().back();
  avNDb = subCollParmsPtr->at(nParms()).data().back();

  // Done.
  return true;
}

//--------------------------------------------------------------------------

// Update the parameters to the interpolated value at the given eCM.

bool SubCollisionModel::setKinematics(double eCMIn) {
  if ( eCMIn > 1.001*eMax || eCMIn < 0.999*eMin ) return false;
  eSave = eCMIn;
  eCMIn = max(eMin, min(eCMIn, eMax));
  if (nParms() > 0) {
    vector<double> parmsNow(subCollParmsPtr->size());
    for (size_t iParm = 0; iParm < parmsNow.size(); ++iParm)
      parmsNow[iParm] = subCollParmsPtr->at(iParm).at(eCMIn);
    avNDb = subCollParmsPtr->at(nParms()).at(eCMIn);
    setParm(parmsNow);
  }
  return true;
}

//--------------------------------------------------------------------------

bool SubCollisionModel::setIDA(int idA) {
  if (nParms() == 0) return true;
  updateSig(idA, idBSave, eSave);
  subCollParmsPtr = &subCollParmsMap[idA];
  idASave = idA;
  if ( sigTarg[0] - sigTarg[6] < 0.001 ) return false;
  return setKinematics(eSave);
}

//--------------------------------------------------------------------------

// Update internally stored cross sections, which in Angantyr should have
// units of femtometer^2.

void SubCollisionModel::updateSig(int idAIn, int idBIn, double eCMIn) {

  if ( sigCmbPtr && eCMIn < eCMlow ) {
    lowEnergyCache.set(sigTargNN, idAIn, idBIn, eCMIn);
    sigTarg = sigTargNN[0];
  }
  else if ( sigCmbPtr ) {
    // Loop over both protons and neutrons.
    for ( int i = 0; i < 4; ++i ) {
      int idA = idAIn;
      int idB = idBIn;
      if ( i%2 == 1 ) {
        if ( projPtr->A() < 2 ) continue;
        idA = ( projPtr->idN() > 0? 2112: -2112 );
      }
      if ( i/2 == 1 ) {
        if ( targPtr->A() < 2 ) continue;
        idB = ( targPtr->idN() > 0? 2112: -2112 );
      }
      double mA = particleDataPtr->m0(idA);
      double mB = particleDataPtr->m0(idB);
      // Total.
      sigTargNN[i][0] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 0)*MB2FMSQ;
      // Non-diffractive.
      sigTargNN[i][1] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 1)*MB2FMSQ;
      // Doubly diffractive.
      sigTargNN[i][2] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 5)*MB2FMSQ;
      // Diffractive (and wounded) projectile.
      sigTargNN[i][3] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 4)*MB2FMSQ +
        sigTargNN[i][1] + sigTargNN[i][2];
      // Diffractive (and wounded) target.
      sigTargNN[i][4] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 3)*MB2FMSQ +
        sigTargNN[i][1] + sigTargNN[i][2];
      // Central diffractive.
      sigTargNN[i][5] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 6)*MB2FMSQ;
      // Elastic.
      sigTargNN[i][6] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 2)*MB2FMSQ;
      // b-slope not used for low energy.
      sigTargNN[i][7] = 0.0;
      // Low energy excitation.
      sigTargNN[i][8] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 7)*MB2FMSQ;
      // Low energy annihilation.
      sigTargNN[i][9] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 8)*MB2FMSQ;
      // Low energy resonance.
      sigTargNN[i][10] =
        sigCmbPtr->sigmaPartial(idA, idB, eCMIn, mA, mB, 9)*MB2FMSQ;
    }
    sigTarg = sigTargNN[0];
    sigErr[7] = 0.0;

  }
  else {

    sigTotPtr->calc(idAIn, idBIn, eCMIn);
    sigTarg[0] = sigTotPtr->sigmaTot()*MB2FMSQ;
    sigTarg[1] = sigTotPtr->sigmaND()*MB2FMSQ;
    sigTarg[2] = sigTotPtr->sigmaXX()*MB2FMSQ;
    sigTarg[3] = sigTotPtr->sigmaAX()*MB2FMSQ + sigTarg[1] + sigTarg[2];
    sigTarg[4] = sigTotPtr->sigmaXB()*MB2FMSQ + sigTarg[1] + sigTarg[2];
    sigTarg[5] = sigTotPtr->sigmaAXB()*MB2FMSQ;
    sigTarg[6] = sigTotPtr->sigmaEl()*MB2FMSQ;
    sigTarg[7] = sigTotPtr->bSlopeEl();
    sigTarg[8] = sigTarg[9] = sigTarg[10] = 0.0;
  }

  // preliminarily set average ND impact parameter as if black disk.
  avNDb = settingsPtr->parm("HeavyIon:SigFitDefAvNDb");
  if ( avNDb <= 0 ) avNDb = 2.0 * sqrt(sigTarg[1]/M_PI) / 3.0;
  if ( avNDb <= 0 ) avNDb = 2.0 * sqrt(sigTarg[0]/M_PI) / 3.0;

}

//--------------------------------------------------------------------------

// Calculate the Chi^2 for the cross section that model in a subclass
// tries to model.

double SubCollisionModel::Chi2(const SigEst & se, int npar) const {

  double chi2 = 0.0;
  int nval = 0;
  for ( int i = 0, Nval = se.sig.size(); i < Nval; ++i ) {
    if ( sigErr[i] == 0.0 ) continue;
    ++nval;
    chi2 += pow2(se.sig[i] - sigTarg[i])/
      (se.dsig2[i] + pow2(sigTarg[i]*sigErr[i]));
  }
  return chi2/double(max(nval - npar, 1));
}


//--------------------------------------------------------------------------

// Helper function to print out stuff.

static void printFit(string name, double fit, double sig, double sigerr,
                 string unit = "mb    ") {
  cout << " |" << setw(25) << name << ": " << setw(8);
  if ( fit >= 100000 )
    cout << "unstable";
  else
    cout << fit;
  cout << (sigerr > 0.0? " *(": "  (")
       << setw(6) << sig;
  if ( sigerr > 0.0 )
    cout << " +- " << setw(2) << int(100.0*sigerr)  << "%";
  else
    cout << "       ";
  cout << ") " << unit << "          | " << endl;
}


//--------------------------------------------------------------------------

// A simple genetic algorithm for fitting the parameters in a subclass
// to reproduce desired cross sections.

bool SubCollisionModel::evolve(int nGenerations, double eCM, int idANow) {

  static int loop = 0;

  if (nParms() == 0)
    return true;
  if (nGenerations <= 0)
    return true;

  if ( fitPrint ) {
    ostringstream os;
    os << "Fitting parameters for " << idANow << " on "  << idBSave
       << " @ " << setprecision(1) << fixed << eCM  << " GeV";
    cout << " |                                      "
         << "                               | \n"
         << " |   " << left << setw(66) << os.str() << "| \n |   ";
    flush(cout);
  }

  // We're going to use a home-made genetic algorithm. We start by
  // creating a population of random parameter points.
  typedef vector<double> Parms;
  Parms minp = minParm();
  Parms maxp = maxParm();
  Parms defp = getParm();
  int dim = nParms();

  // Population of parameter sets. The most accurate sets will propagate
  // to the next generation.
  vector<Parms> pop(NPop, Parms(dim));
  for ( int j = 0; j < dim; ++j )
    pop[0][j] = clamp(defp[j], minp[j], maxp[j]);
  for ( int i = 1; i < NPop; ++i )
    for ( int j = 0; j < dim; ++j )
      pop[i][j] = minp[j] + rndmPtr->flat()*(maxp[j] - minp[j]);

  // Now we evolve our population for a number of generations.
  for ( int iGen = 0; iGen < nGenerations; ++iGen ) {

    // Calculate Chi2 for each parameter set and order them.
    multimap<double, Parms> chi2map;
    double chi2max = 0.0;
    for ( int i = 0; i < NPop; ++i ) {
      setParm(pop[i]);
      double chi2 = Chi2(getSig(), dim);
      chi2map.insert(make_pair(chi2, pop[i]));
      chi2max = max(chi2max, chi2);
    }

    if (fitPrint) {
      if ( iGen >= nGenerations - 20 ) cout << ".";
      flush(cout);
    }

    // Keep the best one, and move the other closer to a better one or
    // kill them if they are too bad.
    multimap<double, Parms>::iterator it = chi2map.begin();
    pop[0] = it->second;
    for ( int i = 1; i < NPop; ++i ) {
      ++it;
      double chi2Now = it->first;
      const Parms& parmsNow = it->second;
      pop[i] = it->second;
      if ( chi2Now > rndmPtr->flat()*chi2max ) {
        // Kill this individual and create a new one.
        for ( int j = 0; j < dim; ++j )
          pop[i][j] = minp[j] + rndmPtr->flat()*(maxp[j] - minp[j]);
      } else {
        // Pick one of the better parameter sets and move this closer.
        int ii = int(rndmPtr->flat()*i);
        for ( int j = 0; j < dim; ++j ) {
          double d = pop[ii][j] - parmsNow[j];
          double pl = clamp(parmsNow[j] - sigFuzz*d, minp[j], maxp[j]);
          double pu = clamp(parmsNow[j] + (1.0 + sigFuzz)*d, minp[j], maxp[j]);
          pop[i][j] = pl + rndmPtr->flat()*(pu - pl);
        }
      }
    }
  }

  // Update resulting parameter set.
  setParm(pop[0]);
  SigEst se = getSig();
  double chi2 = Chi2(se, dim);

  // If the user has deemed the Chi2 too high, continue fitting.
  if ( settingsPtr->parm("HeavyIon:SigFitMaxChi2") > 0.0 &&
       settingsPtr->parm("HeavyIon:SigFitMaxChi2") < chi2 &&
       loop <  settingsPtr->mode("HeavyIon:SigFitMaxChi2Max") ) {
    cout << " Chi2 not converging, continuing...           | \n";
    ++loop;
    bool ret = evolve(nGenerations, eCM, idANow);
    --loop;
    return ret;
  }

  // Output information.
  avNDb = se.avNDb;
  if ( fitPrint ) {
    for ( int i = nGenerations; i < 20; ++i ) cout << " ";
    cout << "                                              | \n";
    cout << " |                                      "
         << "                               | "
         << endl;
    if ( nGenerations > 0 ) {
      cout << " |     Resulting parameters:         "
           << "                                  | "
           << endl;
      for (int iParm = 0; iParm < this->nParms(); ++iParm) {
        cout << right
             << " |" << setw(25) << "[" + to_string(iParm) + "]: "
             << setprecision(2) << setw(7) << pop[0][iParm]
             << setw(39) << "| " << endl;
      }
      cout << " |                                      "
           << "                               | "
           << endl;
    }
    cout << " |     Resulting cross sections        (target value) "
         << "                 | "
         << endl;
    printFit("Total", se.sig[0]*FMSQ2MB,
             sigTarg[0]*FMSQ2MB, sigErr[0]);
    printFit("non-diffractive", se.sig[1]*FMSQ2MB,
             sigTarg[1]*FMSQ2MB, sigErr[1]);
    printFit("XX diffractive", se.sig[2]*FMSQ2MB,
             sigTarg[2]*FMSQ2MB, sigErr[2]);
    printFit("wounded target (B)", se.sig[3]*FMSQ2MB,
             sigTarg[3]*FMSQ2MB, sigErr[3]);
    printFit("wounded projectile (A)", se.sig[4]*FMSQ2MB,
             sigTarg[4]*FMSQ2MB, sigErr[4]);
    printFit("AXB diffractive", se.sig[5]*FMSQ2MB,
             sigTarg[5]*FMSQ2MB, sigErr[5]);
    printFit("elastic", se.sig[6]*FMSQ2MB,
             sigTarg[6]*FMSQ2MB, sigErr[6]);
    printFit("elastic b-slope", se.sig[7], sigTarg[7], sigErr[7], "GeV^-2");
    cout << " |                                   "
         << "                                  | "
         << endl;
    cout << " |" << setw(25) << "Chi2/Ndf" << ": ";
    cout << fixed << setprecision(2);
    cout << setw(8) << chi2 << "                                  | \n";

    cout << " |                                      "
         << "                               | "
         << endl;
  }

  // Done.
  return true;

}

//==========================================================================

// The BlackSubCollisionModel uses fixed size, black-disk
// nucleon-nucleon cross section, equal to the total inelastic pp cross
// section. Everything else is elastic -- Diffraction not included.

//--------------------------------------------------------------------------

SubCollisionSet BlackSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;
  double favNDb = avNDb*impactFudge;

  // Go through all pairs of nucleons
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;
      if ( b < sqrt((sigTot() - sigEl())/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::ABS));
      }
      else {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0,
            SubCollision::ELASTIC));
      }
    }

  return SubCollisionSet(ret, 0.5);
}

//==========================================================================

// The NaiveSubCollisionModel uses a fixed size, black-disk-like
// nucleon-nucleon cross section where. Central collisions will always
// be absorptive, less central will be doubly diffractive, more
// peripheral will be single diffractive and the most peripheral will
// be elastic.

//--------------------------------------------------------------------------

SubCollisionSet NaiveSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;
  double favNDb = avNDb*impactFudge;

  double S = 1.0;

  int sigid = 0;
  for (Nucleon& p : proj) {
    for (Nucleon& t : targ) {

      if ( sigid ) swap(sigTarg, sigTargNN[sigid]);
      sigid = 0;
      if ( projPtr->A() > 1 && abs(p.id()) == 2112 ) sigid += 1;
      if ( targPtr->A() > 1 && abs(t.id()) == 2112 ) sigid += 2;
       if ( sigid ) swap(sigTarg, sigTargNN[sigid]);

      double b = (p.bPos() - t.bPos()).pT();
      if ( b > sqrt(sigTot()/M_PI) ) continue;

      S *= 0.5;

      double currXS = sigND();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::ABS));
        continue;
      }
      currXS += sigLAnn();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::LANN));
        continue;
      }
      currXS += sigLRes();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::LRES));
        continue;
      }
      currXS += sigDDE();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::DDE));
        continue;
      }
      currXS += sigSDE();
      if ( b < sqrt(currXS/M_PI) ) {
        if ( sigSDEP() > rndmPtr->flat()*sigSDE() ) {
          ret.insert(SubCollision(p, t, b, b/favNDb, -1.0,
              SubCollision::SDEP));
        } else {
          ret.insert(SubCollision(p, t, b, b/favNDb, -1.0,
              SubCollision::SDET));
        }
        continue;
      }
      currXS += sigLExc();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::LEXC));
        continue;
      }
      currXS += sigCDE();
      if ( b < sqrt(currXS/M_PI) ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0, SubCollision::CDE));
        continue;
      }
      if ( elasticMode ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, -1.0,
            SubCollision::ELASTIC));
      }
    }
  }
  if ( sigid ) swap(sigTarg, sigTargNN[sigid]);

  return SubCollisionSet(ret, 1.0 - S);
}

//==========================================================================

// DoubleStrikman uses a fluctuating and semi-transparent disk for
// each Nucleon in a sub-collision resulting in a fluctuating
// interaction probability. To assess the fluctuation each Nucleon has
// two random states in each collision, one main state and one helper
// state to assess the fluctuations.

//--------------------------------------------------------------------------

// Virtual init method.
bool FluctuatingSubCollisionModel::init(int idAIn, int idBIn, double eCMIn) {
  return SubCollisionModel::init(idAIn, idBIn, eCMIn);
}

//--------------------------------------------------------------------------

// Helper functions to get the correct average elastic and wounded
// cross sections for fluctuating models.

static void shuffle(double PND1, double PND2, double & PW1, double & PW2) {
  if ( PND1 > PW1 ) {
    PW2 += PW1 - PND1;
    PW1 = PND1;
    return;
  }
  if ( PND2 > PW2 ) {
    PW1 += PW2 - PND2;
    PW2 = PND2;
    return;
  }
}

static void shuffle(double & PEL11, double P11,
                    double P12, double P21, double P22) {
  double PEL12 = PEL11, PEL21 = PEL11, PEL22 = PEL11;
  std::array<std::pair<double, double *>, 4> arr{
      {{P11, &PEL11}, {P12, &PEL12}, {P21, &PEL21}, {P22, &PEL22}}};

  // Sorted in place instead
  std::sort(std::begin(arr), std::end(arr),
            [](const std::pair<double, double *> &a,
               const std::pair<double, double *> &b) {
              return a.first < b.first; // Sort by P values
            });

  for (int i = 0; i < 3; ++i) {
    if (*(arr[i].second) > arr[i].first) {
      *(arr[i + 1].second) += *(arr[i].second) - arr[i].first;
      *(arr[i].second) = arr[i].first;
    }
  }
}

static double pnw(double PWp, double PWt, double PND) {
  return ( 1.0 - PWp <= 0.0 || 1.0 - PWt <= 0.0 )?
    0.0: (1.0 - PWp)*(1.0 - PWt)/(1.0 - PND);
}

static double el(double s1, double s2, double u1, double u2) {
  return s1*u2 > s2*u1? s2*u1: s1*u2;
}

//--------------------------------------------------------------------------

// Numerically estimate the NN semi-inclusive cross sections
// corresponding to the current parameter setting. The radii are
// sampled each iteration, while the integral over impact parametr is
// done analytaically

SubCollisionModel::SigEst FluctuatingSubCollisionModel::getSig() const {

  // FPE prevention.
  const double HUGEVAL = 1.0e100;

  // The random sampling of radii is a bit time-consuming. If we
  // generate a fair amount of them from start and then randomly pick
  // from them for each iteration, we gain a factor of more than 2 in
  // speed.
  int nSample = NInt/4;
  vector<double> projSample(nSample), targSample(nSample);
  for ( int i = 0; i < nSample; ++i ) {
    projSample[i] = pickRadiusProj();
    targSample[i] = pickRadiusTarg();
  }

  SigEst s;
  for ( int n = 0; n < NInt; ++n ) {

    // First we pick 2x2 statistically independent radii combinations.
    double rp1 = projSample[int(nSample*rndmPtr->flat())];
    double rp2 = projSample[int(nSample*rndmPtr->flat())];
    double rt1 = targSample[int(nSample*rndmPtr->flat())];
    double rt2 = targSample[int(nSample*rndmPtr->flat())];
    double s11 = pow2(rp1 + rt1)*M_PI;
    double s12 = pow2(rp1 + rt2)*M_PI;
    double s21 = pow2(rp2 + rt1)*M_PI;
    double s22 = pow2(rp2 + rt2)*M_PI;

    // Calculate the total cross section.
    double stot = (s11 + s12 + s21 + s22)/4.0;
    s.sig[0] += stot;
    s.dsig2[0] += pow2(stot);

    // Calculate the non-diffractive cross section and collect
    // information about corresponding overlap and average impact
    // parameter.
    double u11 = opacity(s11)/2.0;
    double u12 = opacity(s12)/2.0;
    double u21 = opacity(s21)/2.0;
    double u22 = opacity(s22)/2.0;

    if ( s11 < u11*HUGEVAL && s12 < u12*HUGEVAL &&
         s21 < u21*HUGEVAL && s22 < u22*HUGEVAL ) {
      double avb = sqrt(2.0/M_PI)*(s11*sqrt(s11/(2.0*u11))*(1.0 - u11) +
                                   s12*sqrt(s12/(2.0*u12))*(1.0 - u12) +
                                   s21*sqrt(s21/(2.0*u21))*(1.0 - u21) +
                                   s22*sqrt(s22/(2.0*u22))*(1.0 - u22))/12.0;
      s.avNDb += avb;
      s.davNDb2 += pow2(avb);

      double avOF = (0.5*log1mT0(s11)*s11/u11 +
                     0.5*log1mT0(s12)*s12/u12 +
                     0.5*log1mT0(s21)*s21/u21 +
                     0.5*log1mT0(s22)*s22/u22)/4.0;
      s.avNOF +=avOF;
      s.davNOF2 +=pow2(avOF);
    }

    double snd = (s11 - s11*u11 + s12 - s12*u12 +
                  s21 - s21*u21 + s22 - s22*u22)/4.0;
    s.sig[1] += snd;
    s.dsig2[1] += pow2(snd);

    // Calculate the elstic cross section.
    double sel = (el(s11, s22, u11, u22) + el(s12, s21, u12, u21))/2.0;
    s.sig[6] += sel;
    s.dsig2[6] += pow2(sel);

    // Calculate the cross sections for wounded projectile and target.
    double swt = stot - (el(s11, s12, u11, u12) + el(s21, s22, u21, u22))/2.0;
    double swp = stot - (el(s11, s21, u11, u21) + el(s12, s22, u12, u22))/2.0;
    s.sig[4] += swp;
    s.dsig2[4] += pow2(swp);
    s.sig[3] += swt;
    s.dsig2[3] += pow2(swt);

    // Calculate the doubly diffracted cross section.
    s.sig[2] += swt + swp - snd  + sel - stot;
    s.dsig2[2] += pow2(swt + swp - snd  + sel - stot);

    // Calculate the elastic b-slope.
    s.sig[5] += s11;
    s.dsig2[5] += pow2(s11);

    if ( s11*s11 > u11*HUGEVAL ) continue;
    s.sig[7] += pow2(s11)/u11;
    s.dsig2[7] += pow2(pow2(s11)/u11);

  }

  // Normalise everything.
  s.sig[0] /= double(NInt);
  s.dsig2[0] = (s.dsig2[0]/double(NInt) - pow2(s.sig[0]))/double(NInt);

  s.sig[1] /= double(NInt);
  s.dsig2[1] = (s.dsig2[1]/double(NInt) - pow2(s.sig[1]))/double(NInt);

  s.sig[2] /= double(NInt);
  s.dsig2[2] = (s.dsig2[2]/double(NInt) - pow2(s.sig[2]))/double(NInt);

  s.sig[3] /= double(NInt);
  s.dsig2[3] = (s.dsig2[3]/double(NInt) - pow2(s.sig[3]))/double(NInt);

  s.sig[4] /= double(NInt);
  s.dsig2[4] = (s.dsig2[4]/double(NInt) - pow2(s.sig[4]))/double(NInt);

  s.sig[6] /= double(NInt);
  s.dsig2[6] = (s.dsig2[6]/double(NInt) - pow2(s.sig[6]))/double(NInt);

  s.sig[5] /= double(NInt);
  s.dsig2[5] /= double(NInt);

  s.sig[7] /= double(NInt);
  s.dsig2[7] /= double(NInt);

  // Protect from FPEs.
  if ( s.sig[5] > 0.0 && s.sig[7] < s.sig[5]*HUGEVAL ) {
    double bS = (s.sig[7]/s.sig[5])/(16.0*M_PI*pow2(0.19732697));
    double b2S = pow2(bS)*(s.dsig2[7]/pow2(s.sig[7]) - 1.0 +
                           s.dsig2[5]/pow2(s.sig[5]) - 1.0)/double(NInt);
    s.sig[7] = bS;
    s.dsig2[7] = b2S;
  } else {
    s.sig[7] = 0.0;
    s.dsig2[7] = 0.0;
  }
  // We don't really know how to calculate central diffraction yet.
  s.sig[5] = 0.0;
  s.dsig2[5] = 0.0;

  s.avNDb /= double(NInt);
  s.davNDb2 = (s.davNDb2/double(NInt) - pow2(s.avNDb))/double(NInt);

  // Protect from FPEs.
  if ( s.sig[1] > 0.0 ) {
    s.avNDb   /= s.sig[1];
    s.davNDb2 /= pow2(s.sig[1]);
    s.avNOF /= double(NInt);
    s.avNOF /= s.sig[1];
    s.davNOF2 /= double(NInt);
    s.davNOF2 /= pow2(s.sig[1]);
  } else {
    s.avNDb   = 0.0;
    s.davNDb2 = 0.0;
    s.avNOF   = 0.0;
    s.davNOF2  = 0.0;
  }
  return s;

}

//--------------------------------------------------------------------------

// Helper function Given 2x2 statistically equivalent elastic
// amplitudes, shuffle probabilities between them so that the
// different probabilities for inelastic scattering are above zero and
// below unity for all four cases.
//
// Return a vector with the probabilities for non-diffractive, double
// diffractive excitation, single projectile excitation, single target
// excitation, and elastic, for the main amplitude T11. Also trurn
// theprobability of inelastic scatterning for the four statistically
// equivalent amplitudes.

vector<double> FluctuatingSubCollisionModel::
getCollTypeProbs(const vector<double> & T) const {

  if ( T[0] + T[1] + T[2] + T[3] == 0.0 )
    return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  vector<double> PND(5), PNND(5), PNIN(5), PTD(5), PDp(5), PDt(5),
    PDD(5), PEl(5);

  auto op  = [](int j) { return j^2; };
  auto ot  = [](int j) { return j^1; };
  auto opt  = [](int j) { return j^3; };

  // First calculate the differenct dSigma/d2b for the different
  // amplitudes.
  for ( int i = 0; i < 4; ++i ) {
    PND[4] += (PND[i] = 1.0 - pow2(1.0 - T[i]));
    PNND[4] += (PNND[i] = 1.0 - PND[i]);
    PTD[4] += (PTD[i] = pow2(T[i]) - T[i]*T[opt(i)]);
    PDp[4] += (PDp[i] = T[i]*T[ot(i)] - T[i]*T[opt(i)]);
    PDt[4] += (PDt[i] = T[i]*T[op(i)] - T[i]*T[opt(i)]);
    PEl[4] += (PEl[i] = T[i]*T[opt(i)]);
    PDD[4] += (PDD[i] = PTD[i] - PDp[i] - PDt[i]);
  }

  // Then spread out the summed probabilities for diffractive
  // scatterings in proportion to the available probabilities.
  for (int i = 0; i < 4; ++i) {
    PTD[i] = PTD[4]*PNND[i]/PNND[4];
    PDp[i] = PDp[4]*PNND[i]/PNND[4];
    PDt[i] = PDt[4]*PNND[i]/PNND[4];
    PDD[i] = PDD[4]*PNND[i]/PNND[4];
    PNIN[4] += (PNIN[i] = PNND[i] - PTD[i]);
  }

  // Finally see if there is room for elastic scatterings.
  for ( int i = 0; i < 4; ++i )
    if ( PNIN[4] > 0.0 ) PEl[i] = PEl[4]*PNIN[i]/PNIN[4];

  // Now everything (except the elastic probability) should be fine.
  return { PND[0], PDD[0], PDp[0], PDt[0], PEl[0],
    PND[0] + PTD[0], PND[1] + PTD[1], PND[2] + PTD[2], PND[3] + PTD[3] };

}

//--------------------------------------------------------------------------

// Generate radii for all nucleons.

void FluctuatingSubCollisionModel::
generateNucleonStates(Nucleus& proj, Nucleus& targ) {
  // Assign two states to each nucleon.
  for (Nucleon& p : proj) {
    p.state({ pickRadiusProj() });
    p.addAltState({ pickRadiusProj() });
  }
  for (Nucleon& t : targ) {
    t.state({ pickRadiusTarg() });
    t.addAltState({ pickRadiusTarg() });
  }
}

//--------------------------------------------------------------------------

// Main function returning the possible sub-collisions.

SubCollisionSet FluctuatingSubCollisionModel::
getCollisions(Nucleus& proj, Nucleus& targ) {

  if ( opacityMode > 1 ) return getCollisionsNew(proj, targ);

  multiset<SubCollision> ret;
  double favNDb = avNDb*impactFudge;


  // The factorising S-matrix.
  double SS11 = 1.0, SS12 = 1.0, SS21 = 1.0, SS22 = 1.0;

  // Go through all pairs of nucleons
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();
      double olapp = getOverlap(b, t.state()[0],p.state()[0])/avNDolap;
      double T11 = Tpt(p.state(), t.state(), b);
      double T12 = Tpt(p.state(), t.altState(), b);
      double T21 = Tpt(p.altState(), t.state(), b);
      double T22 = Tpt(p.altState(), t.altState(), b);

      double S11 = 1.0 - T11;
      double S12 = 1.0 - T12;
      double S21 = 1.0 - T21;
      double S22 = 1.0 - T22;

      SS11 *= S11;
      SS12 *= S12;
      SS21 *= S21;
      SS22 *= S22;
      double PND11 = 1.0 - pow2(S11);
      // First and most important, check if this is an absorptive
      // scattering.
      if ( PND11 > rndmPtr->flat() ) {
        auto abstype = SubCollision::ABS;
        if ( sigLAnn() + sigLRes() > 0 ) {
          double totabs = (sigND() + sigLAnn() + sigLRes())*rndmPtr->flat();
          if ( totabs > sigND() + sigLAnn() )
            abstype = SubCollision::LRES;
          else if ( totabs > sigND() )
            abstype = SubCollision::LANN;
        }
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp, abstype));
        continue;
      }

      // Now set up calculation for probability of diffractively
      // wounded nucleons.
      double PND21 = 1.0 - pow2(S21);
      double PWp11 = 1.0 - S11*S21;
      double PWp21 = 1.0 - S11*S21;
      shuffle(PND11, PND21, PWp11, PWp21);
      double PND12 = 1.0 - pow2(S12);
      double PWt11 = 1.0 - S11*S12;
      double PWt12 = 1.0 - S11*S12;
      shuffle(PND11, PND12, PWt11, PWt12);

      bool wt = ( PWt11 - PND11 > (1.0 - PND11)*rndmPtr->flat() );
      bool wp = ( PWp11 - PND11 > (1.0 - PND11)*rndmPtr->flat() );
      if ( wt && wp ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp, SubCollision::DDE));
        continue;
      }
      if ( wt ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp, SubCollision::SDET));
        continue;
      }
      if ( wp ) {
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp, SubCollision::SDEP));
        continue;
      }

      // Finally set up calculation for elastic scattering. This can
      // never be exact, but let's do as well as we can.
      if ( elasticMode < 1 ) continue;
      double PND22 = 1.0 - pow2(S22);
      double PWp12 = 1.0 - S12*S22;
      double PWp22 = 1.0 - S12*S22;
      shuffle(PND12, PND22, PWp12, PWp22);
      double PWt21 = 1.0 - S21*S22;
      double PWt22 = 1.0 - S21*S22;
      shuffle(PND21, PND22, PWt21, PWt22);

      double PNW11 = pnw(PWp11, PWt11, PND11);
      double PNW12 = pnw(PWp12, PWt12, PND12);
      double PNW21 = pnw(PWp21, PWt21, PND21);
      double PNW22 = pnw(PWp22, PWt22, PND22);

      double PEL = (T12*T21 + T11*T22)/2.0;
      shuffle(PEL, PNW11, PNW12, PNW21, PNW22);
      if ( PEL > PNW11*rndmPtr->flat() ) {
        if ( sigCDE() > rndmPtr->flat()*(sigCDE() + sigEl()) )
          ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
              SubCollision::CDE));
        else
          ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
              SubCollision::ELASTIC));
      }
    }

  double T11 = 1.0 - SS11;
  double T12 = 1.0 - SS12;
  double T21 = 1.0 - SS21;
  double T22 = 1.0 - SS22;

  return SubCollisionSet(ret, T11, T12, T21, T22);

}

//--------------------------------------------------------------------------

// Main function returning the possible sub-collisions (new version
// for better reproduction of inelastic cross sections).

SubCollisionSet FluctuatingSubCollisionModel::
getCollisionsNew(Nucleus& proj, Nucleus& targ) {

  multiset<SubCollision> ret;
  double favNDb = avNDb*impactFudge;

  // The factorising S-matrix.
  double SS11 = 1.0, SS12 = 1.0, SS21 = 1.0, SS22 = 1.0;

  // Probability of no inelastic scattering.
  double PNI11 = 1.0, PNI12 = 1.0, PNI21 = 1.0, PNI22 = 1.0;

  // Go through all pairs of nucleons.
  for (Nucleon& p : proj)
    for (Nucleon& t : targ) {
      double b = (p.bPos() - t.bPos()).pT();
      double olapp = getOverlap(b, t.state()[0],p.state()[0])/avNDolap;
      double T11 = Tpt(p.state(), t.state(), b);
      double T12 = Tpt(p.state(), t.altState(), b);
      double T21 = Tpt(p.altState(), t.state(), b);
      double T22 = Tpt(p.altState(), t.altState(), b);
      SS11 *= (1.0 - T11);
      SS12 *= (1.0 - T12);
      SS21 *= (1.0 - T21);
      SS22 *= (1.0 - T22);

      auto P = getCollTypeProbs({T11, T12, T21, T22});

      PNI11 *= 1.0 - P[5];
      PNI12 *= 1.0 - P[6];
      PNI21 *= 1.0 - P[7];
      PNI22 *= 1.0 - P[8];
      double R = rndmPtr->flat();
      double acc = 0.0;
      if ( R < (acc += P[0]) )
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
            SubCollision::ABS));
      else if ( R < (acc += P[1]) )
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
            SubCollision::DDE));
      else if ( R < (acc += P[2]) )
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
            SubCollision::SDEP));
      else if ( R < (acc += P[3]) )
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
            SubCollision::SDET));
      else if ( elasticMode > 0 && R < (acc + P[4]) )
        ret.insert(SubCollision(p, t, b, b/favNDb, olapp,
            SubCollision::ELASTIC));
    }

  double T11 = 1.0 - SS11;
  double T12 = 1.0 - SS12;
  double T21 = 1.0 - SS21;
  double T22 = 1.0 - SS22;

  if ( !ret.empty() || elasticMode > 0 || elasticFudge <= 0.0 ||
       proj.size() + targ.size() == 2)
    return SubCollisionSet(ret, T11, T12, T21, T22);

  // Option to include the probability that the nuclei were wounded by
  // elastic NN scatterings.  Calculate the overall probability that
  // the AA collision was inelastic (summed over the four state
  // combinations.
  double PInelSum = 2.0*(T11 + T12 + T21 + T22 - T11*T22 - T12*T21);

  // Now calculate the probability that there were inelastic NN
  // scatterings, again summed over states.
  double PInelNN = 4.0 -  PNI11 - PNI12 - PNI21 - PNI22;

  // Shuffle the difference in probabilities (note that it may be
  // negative, so we need a fudge factor) between states and get the
  // share of elastic NN causing inelastic AA for our primary state.
  double PNNel = PNI11*(PInelSum - PInelNN)/(PNI11 + PNI12 + PNI21 + PNI22);
  if ( PNNel > 0 && PNNel*elasticFudge > PNI11*rndmPtr->flat() ) {
    int idx = int(proj.size()*targ.size()*rndmPtr->flat());
    Nucleon& p = *(proj.begin() + idx/targ.size());
    Nucleon& t = *(targ.begin() + idx%targ.size());
    double b = (p.bPos() - t.bPos()).pT();
    double olapp = getOverlap(b, t.state()[0],p.state()[0])/avNDolap;
    ret.insert(SubCollision(p, t, b, b/favNDb, olapp, SubCollision::ELASTIC));
  }

  return SubCollisionSet(ret, T11, T12, T21, T22);

}

//==========================================================================

} // end namespace Pythia8
