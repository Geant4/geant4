// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.8 1999-07-21 08:40:09 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <rw/tphdict.h>
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

void
G3MatTable::print(G4int matid){
  G4Material* Mat;
  if (matid != -99999) {
    G4cout << "Matid " << matid << " ";
    Mat = get(matid);
    G4cout << Mat << endl;
  } else {
    RWTPtrHashDictionaryIterator<G4String,G4Material> iter(*_Mat);
    if (_Mat->entries() > 0) {
      for (;iter();){
	G4String a = *iter.key();
	G4int lngth = a.length()-3;
	Mat = iter.value();
	G4cout << "Matid " << a(3,lngth) << " ";
	G4cout << Mat << endl;
      }
    } else {
      cerr << "No entries in dictionary" << endl;
    }
  }
}
    





