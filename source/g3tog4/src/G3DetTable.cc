// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTable.cc,v 1.8 1999-12-15 14:49:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3DetTable.hh"

typedef G4std::map<G4String, G3DetTableEntry*, G4std::less<G4String> >
::iterator DTDiterator;

G4String 
G3DetTable::MakeHash(G4String& set, G4String& det){;
  return set+" "+det;
}

G3DetTable::G3DetTable(){
}

G3DetTable::~G3DetTable(){
  if (DTD.size() > 0) {
    //    G4cout << "Deleting DTD" << G4endl;
    for (DTDiterator i=DTD.begin(); i != DTD.end(); i++) {
      delete (*i).second;
    }
    DTD.clear();
  }
};

G4VSensitiveDetector* 
G3DetTable::GetSD(G4String& set, G4String& det){

  // make hash ID
  const G4String ShashID = MakeHash(set, det);

  // search the map
  DTDiterator i = DTD.find(ShashID);
  G3DetTableEntry* DTE = (*i).second;
  if (DTE != 0) {
    return DTE->GetSD();
  } else {
    return 0;
  }  
}

G4int 
G3DetTable::GetID(G4String& set, G4String& det){

  // make hash ID
  G4String ShashID = MakeHash(set, det);

  // search the Hash Dictionary
  DTDiterator i = DTD.find(ShashID);
  G3DetTableEntry* DTE = (*i).second;
  if (DTE != 0) {
    return DTE->GetID();
  } else {
    return 0;
  }
}

void 
G3DetTable::Put(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D){
  // make hash ID
  G4String ShashID = MakeHash(set, det);
  G3DetTableEntry* DTE = new G3DetTableEntry(set, det, id, D);
  G4cout << "Inserted DTE with id " << ShashID << G4endl;
  DTD[ShashID] = DTE;
}

void
G3DetTable::PrintAll(){
  if (DTD.size()>0){
    G4int count=0;
    G4cout << "Dump of DTD - " << DTD.size() << " entries:" << G4endl;
    for (DTDiterator i=DTD.begin(); i != DTD.end(); i++) {
      count++;
      G3DetTableEntry* DTE = (*i).second;
      G4cout << "DTD entry " << G4std::setw(3) << count << " sensitive detector name: " 
	     << DTE->GetSD()->GetName() << G4endl;
    }
  }
}





