// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.cc,v 1.2 1999-05-06 04:22:47 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4strstreambuf.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G3MedTable.hh"

G3MedTable::G3MedTable(){
  _Med = new 
    RWTPtrHashDictionary<RWCString,G4Material>(RWCString::hash);
}

G3MedTable::~G3MedTable(){
  G4cout << "Destructing G3MedTable." << endl;
  _Med->clearAndDestroy();
}

G4Material*
G3MedTable::get(G4int MedID){
  // create a string from the MedID, retrieve pointer from hash dictionary
  G4String Hash(MakeMedID(MedID));
  G4Material* R = _Med->findValue(&Hash);
  if (R == 0) G4cerr << "Arrggh! No G3gstmed found with key " << Hash << endl;
  return R;
};

void 
G3MedTable::put(G4int MedID, G4Material* MatPT){
  // create a string from the MedID, store in hash dictionary
  G4String* Hash = new G4String(MakeMedID(MedID));
  _Med->insertKeyAndValue(Hash, MatPT);
}

char*
G3MedTable::MakeMedID(G4int MedID){
  char buf[20];
  ostrstream ostr(buf, sizeof buf);
  ostr << "Med" << MedID << ends;
  return buf;
}



