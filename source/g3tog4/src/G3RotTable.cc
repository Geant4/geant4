// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.9 1999-11-11 15:35:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <strstream.h>
#include "globals.hh"
#include "G3toG4RotationMatrix.hh"
#include "G3toG4.hh"
#include "G3RotTable.hh"

G3RotTable::G3RotTable(){
  _Rot = new G4RWTPtrHashDictionary<G4String,G3toG4RotationMatrix>(G4String::hash);
};

G3RotTable::~G3RotTable(){
  _Rot->clear();
  delete _Rot;
  // G4cout << "Deleted G3RotTable..." << endl;
};

G3toG4RotationMatrix*
G3RotTable::get(G4int rotid){
  G4String _ShashID; // static
  HashID(rotid, _ShashID);
  return _Rot->findValue(&_ShashID);
};

void 
G3RotTable::put(G4int rotid, G3toG4RotationMatrix* RotPT){
  G4String* _HashID = new G4String(); // Dynamic
  HashID(rotid, _HashID);
  _Rot->insertKeyAndValue(_HashID, RotPT);
};

void
G3RotTable::HashID(G4int rotid, G4String& _HashID){
  char s[20];
  ostrstream ostr(s, sizeof s);
  ostr << "Rot" << rotid << ends;
  _HashID = s;
};

void 
G3RotTable::HashID(G4int rotid, G4String* _HashID){
  HashID(rotid, *_HashID);
};
