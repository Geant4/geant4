#ifndef _G3ROTTABLE_
#define _G3ROTTABLE_ 1

// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.hh,v 1.2 1999-05-06 04:20:44 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/tphdict.h>
#include "G3toG4RotationMatrix.hh"

class G3RotTable {
private:
  RWTPtrHashDictionary<G4String,G3toG4RotationMatrix>* _Rot;
  char* MakeRotID(G4int id);
public:
  G3RotTable();
  ~G3RotTable();
  G3toG4RotationMatrix* get(G4int rotid);
  void put(G4int rotid, G3toG4RotationMatrix* rotpt);
};

extern G3RotTable G3Rot;
#endif
