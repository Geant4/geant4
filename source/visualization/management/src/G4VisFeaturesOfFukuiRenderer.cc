// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisFeaturesOfFukuiRenderer.cc,v 1.2 1999-01-09 16:31:26 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4VisFeaturesOfFukuiRenderer.hh"

G4String G4VisFeaturesOfFukuiRenderer () {
  return
    "High quality technical renderer."
    "\n    Features:      exact hidden line, hidden surface algorithms."
    "\n                   high (unlimited) resolution."
    "\n                   renders to PostScript for viewing and/or hardcopy."
    "\n                   remote rendering."
    "\n                   off-line rendering."
    "\n                   graphical user interface."
    "\n    Disadvantages: compute intensive, takes time (use a fast graphics"
    "\n                   system, such as OpenGL, to select view, then copy"
    "\n                   to this renderer - /vis~/copy/view, /vis~/set/view).";
}
