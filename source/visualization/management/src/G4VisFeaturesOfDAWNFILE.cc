// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisFeaturesOfDAWNFILE.cc,v 1.2 1999-01-09 16:31:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4VisFeaturesOfDAWNFILE.hh"

G4String G4VisFeaturesOfDAWNFILE () {
  return
    "High quality technical renderer."
    "\n    Features:      exact hidden line, hidden surface algorithms."
    "\n                   high (unlimited) resolution."
    "\n                   renders to PostScript for viewing and/or hardcopy."
    "\n                   remote rendering."
    "\n                   off-line rendering."
    "\n                   graphical user interface."
    "\n                   connection via g4.prim file to Fukui Renderer DAWN etc."
    "\n    Disadvantages: compute intensive, takes time (use a fast graphics"
    "\n                   system, such as OpenGL, to select view, then copy"
    "\n                   to this renderer - /vis~/copy/view, /vis~/set/view).";
}
