// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITree.hh,v 1.1 2001-04-10 15:08:46 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy to standard output as
//   ASCII stream.

#ifndef G4ASCIITREE_HH
#define G4ASCIITREE_HH

#include "G4VTree.hh"

class G4ASCIITree: public G4VTree {
public:
  G4ASCIITree ();
  virtual ~G4ASCIITree ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");
};

#endif
