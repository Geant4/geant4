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
//      File name:     G4Data.hh
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 1 Giugno 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "G4Data.hh"


G4Data::~G4Data(){

  this->clear();
}

G4bool G4Data::operator == (const G4Data& input) const {

  return(this->length() == input.length());

}

G4bool G4Data::operator < (const G4Data& input) const {

  return(this->length() < input.length());

}



