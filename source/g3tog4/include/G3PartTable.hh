// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3PartTable.hh,v 1.4 1999-12-05 17:50:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G3 materials table.
// Maps G3 material indices to their G4 material object counterparts.
// Maintains a linked List of G3 material index/G4 material pointer pairs.

#include "g4std/map"
#include "G4ParticleDefinition.hh"

class G3PartTable {
private:
  G4std::map<G4String, G4ParticleDefinition*, less<G4String> > PTD;
  void HashID(G4int partid, G4String* _HID);
  void HashID(G4int partid, G4String& _HID);
public:
  G3PartTable();
  virtual ~G3PartTable();
  G4ParticleDefinition* Get(G4int partid);
  void Put(G4int partid, G4ParticleDefinition* partpt);
  void PrintAll();
};

extern G3PartTable G3Part;
