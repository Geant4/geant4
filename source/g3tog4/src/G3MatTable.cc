// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.7 1999-07-20 14:16:41 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <strstream.h>
#include "globals.hh"
#include "G3MatTable.hh"
#include "G4Material.hh"


G3MatTable::G3MatTable(){
  _Mat = new RWTPtrHashDictionary<G4String,G4Material>(G4String::hash);
};

G3MatTable::~G3MatTable(){
  _Mat->clear();
  G4cout << "Deleted G3MatTable..." << endl;
  delete _Mat;
};

G4Material*
G3MatTable::get(G4int matid){
  G4String _ShashID; // static
  HashID(matid, _ShashID);
  return _Mat->findValue(&_ShashID);
};

void 
G3MatTable::put(G4int matid, G4Material* MatPT){
  G4String* _HashID = new G4String(); // Dynamic
  HashID(matid, _HashID);
  _Mat->insertKeyAndValue(_HashID, MatPT);
};

void
G3MatTable::HashID(G4int matid, G4String& _HashID){
  char s[20];
  ostrstream ostr(s, sizeof s);
  ostr << "Mat" << matid << ends;
  _HashID = s;
};

void 
G3MatTable::HashID(G4int matid, G4String* _HashID){
  HashID(matid, *_HashID);
};
