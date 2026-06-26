// HIBasics.h is a part of the PYTHIA event generator.
// Copyright (C) 2026 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the EventInfo class.
//
// EventInfo: stores full nucleon-nucleon events with corresponding Info.

#ifndef Pythia8_HIBasics_H
#define Pythia8_HIBasics_H

#include "Pythia8/Pythia.h"

namespace Pythia8 {

// Forward declarations.
class Nucleon;
class SubCollision;

//==========================================================================

// Class for storing Events and Info objects.

class EventInfo {

public:

  // Empty constructor.
  EventInfo() = default;

  // The Event object.
  Event event;

  // The corresponding Info object.
  Info info;

  // The code for the subprocess.
  int code = 0;

  // The ordering variable of this event.
  double ordering = -1.0;
  bool operator<(const EventInfo & ei) const {
    return ordering < ei.ordering;
  }

  // The associated SubCollision object.
  const SubCollision* coll = {};

  // Is the event properly generated?
  bool ok = false;

  // Which projectile and target nucleons are included and where are
  // they placed?
  map<Nucleon*, pair<int,int> > projs, targs;

  // Also map the remnants for thie for the projectile and/or target, if any.
  map<Nucleon*, vector<int> > projRems, targRems;

  // ... and the location of the projectile's and/or target's (quasi)
  // elastically scattered nucleon. (Note that this is also a remnant.)
  map<Nucleon*, int> projEl, targEl;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HIBasics_H
