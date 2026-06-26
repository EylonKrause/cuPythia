// RivetHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Philip Ilten.

// This class implements an interface to Rivet 4 that can be loaded
// via the plugin structure. It can be run with PythiaParallel and merges
// the separate analyses into a single one at the end of the run.

#ifndef Pythia8_RivetHooks_H
#define Pythia8_RivetHooks_H

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"
#include "Pythia8Plugins/HepMC3.h"

// Rivet includes.
#include "Rivet/AnalysisHandler.hh"

// Directory creation for POSIX.
#include <sys/stat.h>

namespace Pythia8 {

//==========================================================================

// Create a directory, with nesting if necessary. This is not compatible
// outside POSIX systems, and should be removed after migration to C++17.

bool makeDir(vector<string> path, mode_t mode = 0777) {

  // Create the structure needed for stat.
  struct stat info;

  // Loop over the directories.
  string pathNow = "";
  bool first = true;
  for (string& dir : path) {

    // Check if the directory exists.
    pathNow += (first ? "" : "/") + dir;
    first = false;
    if (stat(pathNow.c_str(), &info) == 0) {
      if ((info.st_mode & S_IFDIR) == 0) return false;
      else continue;
    }

    // Create the directory and check it exists.
    if (mkdir(pathNow.c_str(), mode) != 0) return false;
  }
  return true;

}

//==========================================================================

// UserHook to run Rivet analyses.

class RivetHooks : public UserHooks {

public:

  // Constructors and destructor.
  RivetHooks() {}
  RivetHooks(Pythia* pythiaPtrIn, Settings*, Logger*) :
    pythiaPtr(pythiaPtrIn) {}
  ~RivetHooks() {if (rivetPtr != nullptr) delete rivetPtr;}

  //--------------------------------------------------------------------------

  // Perform the Rivet analysis.
  void onEndEvent(Status) override {

    // Optionally skip zero-weight events.
    if (flag("Rivet:skipZeroWeights") && pythiaPtr->info.weight() == 0) {
      loggerPtr->INFO_MSG("skipping zero-weight event");
      return;
    }

    // Convert the event to HepMC.
    hepmc.fillNextEvent(*pythiaPtr);

    // Create the Rivet analysis handler if it does not exist.
    // The analysis handler constructor is not thread safe, so the
    // global mutex must be locked.
    if (rivetPtr == nullptr) {
      mutexPtr->lock();

      // Create the handler.
      rivetPtr = new Rivet::AnalysisHandler();

      // Set whether beams should be checked.
      rivetPtr->setCheckBeams(flag("Rivet:checkBeams"));

      // Set file dumping if requested.
      if (mode("Rivet:dumpPeriod") > 0) {
        string dumpName = word("Rivet:dumpName");
        if (dumpName == "") dumpName = word("Rivet:fileName");
        rivetPtr->setFinalizePeriod(dumpName, mode("Rivet:dumpPeriod"));
      }

      // Preload requested data.
      for (string& preload : wvec("Rivet:preloads"))
        rivetPtr->readData(preload);

      // Set the analyses.
      for (string& analysis : wvec("Rivet:analyses"))
        rivetPtr->addAnalysis(analysis);

      // Initialize Rivet and unlock the global mutex.
      rivetPtr->init(hepmc.event());
      mutexPtr->unlock();
    }

    // Run the analysis.
    rivetPtr->analyze(hepmc.event());

  }

  //--------------------------------------------------------------------------

  // Write the Rivet analyses.
  void onStat() override {

    // Create the directory structure if needed.
    if (rivetPtr == nullptr) return;
    string filename = word("Rivet:fileName");
    vector<string> path = splitString(filename, "/");
    if (path.size() > 1) {
      mutexPtr->lock();
      makeDir(vector<string>(path.begin(), path.end() - 1));
      mutexPtr->unlock();
    }

    // Write the Rivet output.
    rivetPtr->finalize();
    rivetPtr->writeData(filename);

  }

  //--------------------------------------------------------------------------

  // Merge the Rivet analysis handlers.
  void onStat(vector<PhysicsBase*> hookPtrs, Pythia*) override {

    // Get the main Rivet hook.
    Rivet::AnalysisHandler* rivetMain = rivetPtr;
    if (rivetMain == nullptr) {
      loggerPtr->ERROR_MSG("could not retrieve first RivetHooks");
      return;
    }

    // Get the Rivet pointer for each thread.
    int iPtr = 0;
    for (PhysicsBase* hookPtr : hookPtrs) {
      if (hookPtr == this) continue;
      RivetHooks* hookNow = dynamic_cast<RivetHooks*>(hookPtr);
      if (hookNow == nullptr) {
        loggerPtr->ERROR_MSG("could not retrieve RivetHooks for thread ",
          toString(iPtr));
        return;
      } else ++iPtr;

      // Get the current anaysis handler.
      Rivet::AnalysisHandler* rivetNow = hookNow->rivetPtr;
      if (rivetNow == nullptr) {
        loggerPtr->ERROR_MSG("could not retrieve AnalysisHandler for thread ",
          toString(iPtr));
        return;
      }

      // Merge the analysis handler with the main anaysis handler.
      // This should be done with AnalysisHandler::serializeContent and
      // AnalysisHandler::deserializeContent. However, with the most recent
      // Rivet versions, this appears to cause issues in a threaded
      // environment, and so for now we use use AnalysisHandler::merge.
      try {
        rivetMain->merge(*rivetNow);
      } catch (const Rivet::UserError err) {
        loggerPtr->ERROR_MSG("failed to merge analyses in thread ",
          toString(iPtr));
        return;
      }
      rivetMain->merge(*rivetNow);
    }

    // Write the data with the main Rivet hook.
    onStat();

  }

  //--------------------------------------------------------------------------

 private:

  Pythia* pythiaPtr{};
  Rivet::AnalysisHandler* rivetPtr{};
  Pythia8ToHepMC hepmc;

};

//--------------------------------------------------------------------------

// Register Rivet settings.

void rivetSettings(Settings *settingsPtr) {
  settingsPtr->addWord("Rivet:fileName", "rivet.yoda");
  settingsPtr->addWVec("Rivet:analyses", {});
  settingsPtr->addWVec("Rivet:preloads", {});
  settingsPtr->addFlag("Rivet:checkBeams", true);
  settingsPtr->addWord("Rivet:dumpName", "");
  settingsPtr->addMode("Rivet:dumpPeriod", -1, true, false, -1, 0);
  settingsPtr->addFlag("Rivet:skipZeroWeights", false);
}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(UserHooks, RivetHooks, true, false, false)
PYTHIA8_PLUGIN_SETTINGS(rivetSettings)
PYTHIA8_PLUGIN_PARALLEL(true);
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_RivetHooks_H
