// HepMC3Hooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Christian T. Preuss.

// This class implements an interface to HepMC3 that can be loaded
// via the plugin structure. It can be run with PythiaParallel.

#ifndef Pythia8_HepMC3Hooks_H
#define Pythia8_HepMC3Hooks_H

// Pythia includes.
#include "Pythia8/Pythia.h"
#include "Pythia8/Plugins.h"
#include "Pythia8Plugins/HepMC3.h"

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

// UserHook to write HepMC3 files.

class HepMC3Hooks : public UserHooks {

public:

  // Constructors and destructor.
  HepMC3Hooks() {}
  HepMC3Hooks(Pythia* pythiaPtrIn, Settings*, Logger*) :
    pythiaPtr(pythiaPtrIn) {}
  ~HepMC3Hooks() {if (hepMCPtr != nullptr) delete hepMCPtr;}

  //--------------------------------------------------------------------------

  // Print event to HepMC file.
  void onEndEvent(Status) override {

    // Create the HepMC converter.
    if (hepMCPtr == nullptr) {

      // Set the filename if running in parallel.
      string filename = word("HepMC:filename");
      int idx = mode("Parallelism:index");
      if (idx >= 0) {
        size_t iSuffix = filename.find(".hepmc");
        if (iSuffix != string::npos)
          filename = filename.substr(0,iSuffix);
        filename = filename + "_" + to_string(idx) + ".hepmc";
      }

      // Create a HepMC converter and directory structure if needed.
      vector<string> path = splitString(filename, "/");
      if (path.size() > 1) {
        mutexPtr->lock();
        makeDir(vector<string>(path.begin(), path.end() - 1));
        mutexPtr->unlock();
      }
      hepMCPtr = new Pythia8ToHepMC(filename);

      // Save some settings.
      hepMCPtr->set_print_inconsistency(flag("HepMC:printInconsistency"));
      hepMCPtr->set_free_parton_warnings(flag("HepMC:freePartonWarnings"));
      hepMCPtr->set_store_pdf(flag("HepMC:storePDF"));
      hepMCPtr->set_store_proc(flag("HepMC:storeProcess"));
    }

    // Convert the event to HepMC.
    hepMCPtr->fillNextEvent(*pythiaPtr);

    // Write event. Currently each thread writes into its own file.
    hepMCPtr->writeEvent();

  }

  //--------------------------------------------------------------------------

  // Finalise.
  void onStat() override {}

  //--------------------------------------------------------------------------

  // Finalise. Currently, nothing is done here, but merging of the HepMC3
  // records could be done here if this is a useful feature.
  void onStat(vector<PhysicsBase*>, Pythia*) override {
    onStat();
  }

private:

  Pythia* pythiaPtr{};
  Pythia8ToHepMC* hepMCPtr{};

};

//--------------------------------------------------------------------------

// Register HepMC settings.

void hepmcSettings(Settings *settingsPtr) {
  settingsPtr->addWord("HepMC:fileName", "events.hepmc");
  settingsPtr->addFlag("HepMC:printInconsistency", "true");
  settingsPtr->addFlag("HepMC:freePartonWarnings", "true");
  settingsPtr->addFlag("HepMC:storePDF", "true");
  settingsPtr->addFlag("HepMC:storeProcess", "true");
}

//--------------------------------------------------------------------------

// Declare the plugin.

PYTHIA8_PLUGIN_CLASS(UserHooks, HepMC3Hooks, true, false, false)
PYTHIA8_PLUGIN_SETTINGS(hepmcSettings)
PYTHIA8_PLUGIN_PARALLEL(true)
PYTHIA8_PLUGIN_VERSIONS(PYTHIA_VERSION_INTEGER)

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_HepMC3Hooks_H
