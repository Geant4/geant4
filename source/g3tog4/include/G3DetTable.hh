// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTable.hh,v 1.2 1999-05-07 04:15:57 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include <rw/tphdict.h>

#include "globals.hh"

class G4VSensitiveDetector;

class DetTableEntry {
private:
  G4String _set;
  G4String _det;
  G4int _id;
  G4VSensitiveDetector* _detpt;

public:
  DetTableEntry(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D);
  ~DetTableEntry();
  G4VSensitiveDetector* getSD();
  G4String getset();
  G4String getdet();
  G4int getid();
};

class G3DetTable {
private:
  RWTPtrHashDictionary<G4String, DetTableEntry>* _Det;
  G4String MakeHash(G4String& set, G4String& det);
  DetTableEntry* _DTE;

public:
  G3DetTable();
  ~G3DetTable();
  G4VSensitiveDetector* getSD(G4String& set, G4String& det);
  G4int GetID(G4String& set, G4String& det);
  void put(G4String& set, G4String& det, G4int id, G4VSensitiveDetector* D);
  
};

extern G3DetTable G3Det;
