// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.2 1999-05-06 04:23:24 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4strstreambuf.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G3toG4RotationMatrix.hh"
#include "G3RotTable.hh"

G3RotTable::G3RotTable(){
  _Rot = new 
    RWTPtrHashDictionary<RWCString,G3toG4RotationMatrix>(RWCString::hash);
}

G3RotTable::~G3RotTable(){
  G4cout << "Destructing G3RotTable." << endl;
  _Rot->clearAndDestroy();
}

G3toG4RotationMatrix*
G3RotTable::get(G4int RotID){
  // create a string from the RotID, retrieve pointer from hash dictionary
  G4String Hash(MakeRotID(RotID));
  G3toG4RotationMatrix* R = _Rot->findValue(&Hash);
  if (R == 0) {
    G4cerr << "Arrggh! No rotation matrix found with key " << Hash << endl;
  }
  return R;
};

void 
G3RotTable::put(G4int RotID, G3toG4RotationMatrix *RotPT){
  // create a string from the RotID, store in hash dictionary
  G4String* Hash = new G4String(MakeRotID(RotID));
  _Rot->insertKeyAndValue(Hash, RotPT);
}

char*
G3RotTable::MakeRotID(G4int RotID){
  char buf[20];
  ostrstream ostr(buf, sizeof buf);
  ostr << "Rot" << RotID << ends;
  return buf;
}
