// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3DetTable.cc,v 1.3 1999-05-07 04:16:11 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3DetTable.hh"
#include "G4VSensitiveDetector.hh"

DetTableEntry::DetTableEntry(G4String& set, G4String& det, G4int id, 
			     G4VSensitiveDetector* D){
  _set = set;
  _det = det;
  _id  = id;
  _detpt = D;
};

DetTableEntry::~DetTableEntry(){;}

G4VSensitiveDetector* 
DetTableEntry::getSD(){
  return _detpt;
}

G4String 
DetTableEntry::getset(){
  return _set;
}

G4String 
DetTableEntry::getdet(){
  return _det;
}

G4int
DetTableEntry::getid(){
  return _id;
};

G4String 
G3DetTable::MakeHash(G4String& set, G4String& det){;
  return set+" "+det;
};

G3DetTable::G3DetTable(){
  _Det = new RWTPtrHashDictionary<G4String,DetTableEntry>(RWCString::hash);
};

G3DetTable::~G3DetTable(){
  _Det->clearAndDestroy();
  delete _Det;
};

G4VSensitiveDetector* 
G3DetTable::getSD(G4String& set, G4String& det){

  // make hash ID
  G4String HashID = MakeHash(set, det);

  // search the Hash Dictionary
  _DTE = _Det->findValue(&HashID);
  if (_DTE != 0) {
    return _DTE->getSD();
  } else {
    G4cerr << "No G3DetTable entry with key '" << HashID << "'" << endl;
    return 0;
  }
};

G4int 
G3DetTable::GetID(G4String& set, G4String& det){

  // make hash ID
  G4String HashID = MakeHash(set, det);

  // search the Hash Dictionary
  _DTE = _Det->findValue(&HashID);
  if (_DTE != 0) {
    return _DTE->getid();
  } else {
    G4cerr << "No G3DetTable entry with key '" << HashID << "'" << endl;
    return 0;
  }
};

void 
G3DetTable::put(G4String& set, G4String& det, G4int id, 
		G4VSensitiveDetector* D){
  G4String* _HashID = new G4String(MakeHash(set, det));
  _DTE = new DetTableEntry(set, det, id, D);
  _Det->insertKeyAndValue(_HashID, _DTE);
};

