// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3PartTable.hh,v 1.7 2000-11-24 09:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class description:
//
// G3 particles table.
// Maps G3 particles indices to their G4 particle object counterparts.

// ----------------------

#ifndef G3PARTTABLE_HH
#define G3PARTTABLE_HH 1
#include "g4std/map"
#include "G4ParticleDefinition.hh"

class G3PartTable
{

public:  // with description

  G3PartTable();
  virtual ~G3PartTable();
  G4ParticleDefinition* Get(G4int partid);
  void Put(G4int partid, G4ParticleDefinition* partpt);
  void PrintAll();

private:

  G4std::map<G4String, G4ParticleDefinition*, G4std::less<G4String> > PTD;
  void HashID(G4int partid, G4String* _HID);
  void HashID(G4int partid, G4String& _HID);
};

extern G3PartTable G3Part;

#endif
