// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4SecondLevel.hh
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 1 Giugno 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "G4SecondLevel.hh"


G4SecondLevel::~G4SecondLevel(){

  this->clearAndDestroy();
}

G4bool G4SecondLevel::operator == (const G4SecondLevel& input) const{

  return( this->entries() == input.entries());

}

G4bool G4SecondLevel::operator < (const G4SecondLevel& input) const{

  return(this->entries() < input.entries());

}


