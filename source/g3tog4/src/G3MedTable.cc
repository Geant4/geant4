// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MedTable.cc,v 1.5 1999-05-28 21:01:41 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G3MedTable.hh"
#include "G4Material.hh"

G3MedTable::G3MedTable(){
  _Med = new 
    RWTPtrOrderedVector<G4Material>;
};

G3MedTable::~G3MedTable(){
  G4cout << "Clearing G3MedTable." << endl;
  _Med->clear();
  G4cout << "Destructing G3MedTable." << endl;
  delete _Med;
};

G4Material*
G3MedTable::get(G4int MedID){
  return (*_Med)[MedID-1];
};

void 
G3MedTable::put(G4int MedID, G4Material* MatPT){
  _Med->insertAt(MedID-1, MatPT);
};
