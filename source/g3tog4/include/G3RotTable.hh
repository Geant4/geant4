#ifndef _G3ROTTABLE_
#define _G3ROTTABLE_ 1

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.hh,v 1.8 1999-11-23 04:27:25 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Maps G3 rotation indices to G3toG4RotationMatrix*

#include "g4rw/tphdict.h"
#include "G3toG4RotationMatrix.hh"

class G3RotTable {
private:
  G4RWTPtrHashDictionary<G4String, G3toG4RotationMatrix>* _Rot;
  void HashID(G4int rotid, G4String* _HID);
  void HashID(G4int rotid, G4String& _HID);
public:
  G3RotTable();
  ~G3RotTable();
  G4RotationMatrix* get(G4int rotid);
  void put(G4int rotid, G3toG4RotationMatrix* rotpt);
};

extern G3RotTable G3Rot;
#endif
