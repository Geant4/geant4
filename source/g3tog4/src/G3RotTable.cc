// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3RotTable.cc,v 1.6 1999-05-28 21:09:06 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G3toG4RotationMatrix.hh"
#include "G3toG4.hh"
#include "G3RotTable.hh"

G3RotTable::G3RotTable(){
  _Rot = new 
    RWTPtrOrderedVector<G3toG4RotationMatrix>;
};

G3RotTable::~G3RotTable(){
  _Rot->clearAndDestroy();
  delete _Rot;
  G4cout << "Deleted G3RotTable..." << endl;
};

G3toG4RotationMatrix*
G3RotTable::Get(G4int RotID){
  return (*_Rot)[RotID];
};

void 
G3RotTable::Put(G4int RotID, G3toG4RotationMatrix *RotPT){
  _Rot->insertAt(RotID, RotPT);
};
