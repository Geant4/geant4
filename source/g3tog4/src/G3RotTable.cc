// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.4 1999-05-06 18:12:29 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G3toG4RotationMatrix.hh"
#include "G3RotTable.hh"

G3RotTable::G3RotTable(){
  _Rot = new 
    RWTPtrOrderedVector<G3toG4RotationMatrix>;
};

G3RotTable::~G3RotTable(){
  G4cout << "Destructing G3RotTable." << endl;
  _Rot->clearAndDestroy();
};

G3toG4RotationMatrix*
G3RotTable::get(G4int RotID){
  return (*_Rot)[RotID-1];
};

void 
G3RotTable::put(G4int RotID, G3toG4RotationMatrix *RotPT){
  _Rot->insertAt(RotID-1, RotPT);
};
