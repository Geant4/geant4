// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisFeaturesOfRay.cc,v 1.2 1999-01-09 16:31:27 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VisFeaturesOfRay.hh"

G4String G4VisFeaturesOfRayX () {
  return
    "    GEANT4 Ray Tracing to a dumb single buffered X Window."
    "\n    Advantages:    uses GEANT4's own tracking algorithms."
    "\n                   photo-realistic."
    "\n    Disadvantages: slow.";
}
