// PythiaPython.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Aartem Havryliuk and Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to load purely Pythonic files into the pythia8 module.

#ifndef awkward_PythiaPython_H
#define awkward_PythiaPython_H

#include <pybind11/eval.h>

namespace Pythia8 {

//==========================================================================

// Load purely Pythonic modules.

void loadExtraPython(pybind11::module scope);

//==========================================================================

} // end namespace std

#endif // awkward_PythiaPython_H
