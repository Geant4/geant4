// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.5 1999-05-28 20:54:05 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G3MatTable.hh"
#include "G4Material.hh"

G3MatTable::G3MatTable(){
  _Mat = new 
    RWTPtrOrderedVector<G4Material>;
};

G3MatTable::~G3MatTable(){
  G4cout << "Clearing G3MatTable list." << endl;
  _Mat->clear();
  G4cout << "Deleting G3MatTable list." << endl;
  delete _Mat;
};

G4Material*
G3MatTable::get(G4int MatID){
  return (*_Mat)[MatID-1];
};

void G3MatTable::put(G4int MatID, G4Material* MatPT){
  _Mat->insertAt(MatID-1, MatPT);
};
