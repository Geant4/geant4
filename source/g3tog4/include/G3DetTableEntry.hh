// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTableEntry.hh,v 1.3 1999-12-09 01:27:42 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G3DetTableEntry class

#ifndef DETTABLEENTRY_HH
#define DETTABLEENTRY_HH 1

#include "g4std/map"
#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class G3DetTableEntry {
private:
  G4String _set;
  G4String _det;
  G4int _id;
  G4VSensitiveDetector* _detpt;

public:
  G3DetTableEntry(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D);
  ~G3DetTableEntry();
  G4VSensitiveDetector* GetSD();
  G4String GetSet();
  G4String GetDet();
  G4int GetID();
};
#endif
