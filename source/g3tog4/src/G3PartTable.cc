// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3PartTable.cc,v 1.7 1999-12-05 17:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/strstream"
#include "g4std/iomanip"
#include "G3PartTable.hh"

typedef G4std::map<G4String, G4ParticleDefinition*, less<G4String> >
::iterator PTDiterator;

G3PartTable::G3PartTable(){
}

G3PartTable::~G3PartTable(){
  if (PTD.size()>0){
    //    G4cout << "Deleting PTD" << endl;
    for (PTDiterator i=PTD.begin(); i != PTD.end(); i++) {
      delete (*i).second;
    }
    PTD.clear();
  }
}

G4ParticleDefinition*
G3PartTable::Get(G4int partid){
  G4String ShashID; // static
  HashID(partid, ShashID);
  PTDiterator i = PTD.find(ShashID);
  return (*i).second;
};

void 
G3PartTable::Put(G4int partid, G4ParticleDefinition *partpt){
  G4String ShashID; // static
  HashID(partid, ShashID);
  PTD[ShashID]=partpt;
};

void
G3PartTable::HashID(G4int partid, G4String& theHashID){
  char s[20];
  ostrstream ostr(s, sizeof s);
  ostr << "Part" << partid << ends;
  theHashID = s;
};

void 
G3PartTable::HashID(G4int partid, G4String* theHashID){
  HashID(partid, *theHashID);
};

void
G3PartTable::PrintAll(){
  if (PTD.size()>0){
    G4int count=0;
    G4cout << "Dump of PTD - " << PTD.size() << " entries: " << endl;
    for (PTDiterator i=PTD.begin(); i != PTD.end(); i++) {
      count++;
      G4ParticleDefinition* aPTD = (*i).second;
      G4cout << "PTD entry " << setw(3) << count << " particle name: " 
	     << aPTD->GetParticleName() << endl;
    }
  }
}





