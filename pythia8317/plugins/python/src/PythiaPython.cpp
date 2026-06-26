// PythiaPython.cpp is a part of the PYTHIA event generator.
// Copyright (C) 2026 Philip Ilten and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Definitions for loading purely Pythonic files into the pythia8 module.

#include <dirent.h>
#include "extra/PythiaPython.h"
#include "Pythia8/PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================

// Load purely Pythonic modules.

void loadExtraPython(pybind11::module scope) {

  // Find the Python files.
  vector<string> pys;
  string path("python");
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(path.c_str())) != nullptr) {
    while ((ent = readdir (dir)) != nullptr) {
      pys.push_back(ent->d_name);
    }
    closedir(dir);
  }
  
  // Execute the pure Pythonic files.
  for (const string& py : pys) {
    ifstream is(path + "/" + py);
    stringstream ss;
    ss << is.rdbuf();
    pybind11::exec(ss.str(), scope.attr("__dict__"));
  }
  
}

//==========================================================================

} // end namespace Pythia8

