// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.cc,v 1.7 1999-07-20 14:16:45 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <strstream.h>
#include "globals.hh"
#include "G3MedTable.hh"
#include "G4Material.hh"

G3MedTable::G3MedTable(){
  _Med = new RWTPtrHashDictionary<G4String,G4Material>(G4String::hash);
};

G3MedTable::~G3MedTable(){
  _Med->clear();
  G4cout << "Deleted G3MedTable..." << endl;
  delete _Med;
};

G4Material*
G3MedTable::get(G4int medid){
  G4String _ShashID; // static
  HashID(medid, _ShashID);
  return _Med->findValue(&_ShashID);
};

void 
G3MedTable::put(G4int medid, G4Material* MatPT){
  G4String* _HashID = new G4String(); // Dynamic
  HashID(medid, _HashID);
  _Med->insertKeyAndValue(_HashID, MatPT);
};

void
G3MedTable::HashID(G4int medid, G4String& _HashID){
  char s[20];
  ostrstream ostr(s, sizeof s);
  ostr << "Med" << medid << ends;
  _HashID = s;
};

void 
G3MedTable::HashID(G4int medid, G4String* _HashID){
  HashID(medid, *_HashID);
};
