// main225.cc is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Philip Ilten <philten@cern.ch>

// Keywords: command file; rivet; parallelism

// Demonstrates how to use the RivetHooks plugin to run Rivet
// analyses, both with a single Pythia instance and PythiaParallel.

#include "Pythia8/PythiaParallel.h"

using namespace Pythia8;
int main() {

  // Example run with Pythia parallel.
  PythiaParallel pythiaP;
  pythiaP.readFile("main225.cmnd", 0);
  if (!pythiaP.init()) return 1;
  pythiaP.run([&](Pythia*){});
  pythiaP.stat();

  // Same run, but with single Pythia instance.
  Pythia pythiaS;
  pythiaS.readFile("main225.cmnd", 1);
  if (!pythiaS.init()) return 1;
  int nEvent = pythiaS.mode("Main:numberOfEvents");
  int iEvent = 0;
  while (iEvent < nEvent)
    if (pythiaS.next()) ++iEvent;
  pythiaS.stat();
  return 0;

}
