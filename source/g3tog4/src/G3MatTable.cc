// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.2 1999-05-06 04:22:38 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ios.hh"
#include "G4strstreambuf.hh"
#include "globals.hh"
#include "G3MatTable.hh"

G3MatTable::G3MatTable(){
  _Mat = new 
    RWTPtrHashDictionary<RWCString,G4Material>(RWCString::hash);
};

G3MatTable::~G3MatTable(){
  G4cout << "Destructing G3MatTable." << endl;
  _Mat->clearAndDestroy();
  delete _Mat;
};

G4Material*
G3MatTable::get(G4int MatID){
  // create a string from the MatID, retrieve pointer from hash dictionary
  G4String Hash(MakeMatID(MatID));
  G4Material* R = _Mat->findValue(&Hash);
  if (R == 0) G4cerr << "Arrggh! No G4Material found with key " << Hash << endl;
  return R;
};

void G3MatTable::put(G4int MatID, G4Material* MatPT){
  // create a string from the MatID, store in hash dictionary
  G4String* Hash = new G4String(MakeMatID(MatID)); 
  _Mat->insertKeyAndValue(Hash, MatPT);
};

char*
G3MatTable::MakeMatID(G4int MatID){
  char buf[20];
  ostrstream ostr(buf, sizeof buf);
  ostr << "Mat" << MatID << ends;
  return buf;
}
