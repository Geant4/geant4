// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3MatTable.cc,v 1.3 1999-05-06 17:46:53 lockman Exp $
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
  G4cout << "Destructing G3MatTable." << endl;
  _Mat->clearAndDestroy();
  delete _Mat;
};

G4Material*
G3MatTable::get(G4int MatID){
  return (*_Mat)[MatID];
};

void G3MatTable::put(G4int MatID, G4Material* MatPT){
  _Mat->insertAt(MatID, MatPT);
};
