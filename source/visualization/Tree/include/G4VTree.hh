// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VTree.hh,v 1.1 2001-04-10 15:08:48 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#ifndef G4VTREE_HH
#define G4VTREE_HH

#include "G4VGraphicsSystem.hh"

class G4VTree: public G4VGraphicsSystem {
public:
  G4VTree (const G4String& name,
	   const G4String& nickname,
	   const G4String& description,
	   Functionality f);
  virtual ~G4VTree ();
};

#endif
