// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisFeaturesOfDAWNFILE.cc,v 1.4 2001-02-05 02:34:29 johna Exp $
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
    "\n                   view parameters to this renderer:"
    "\n                   \"/vis/viewer/set/all <OpenGL-view>\".";
}
