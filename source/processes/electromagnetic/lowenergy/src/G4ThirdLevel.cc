//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4ThirdLevel.hh
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 1 Giugno 1999
//
//      Modifications: 24.04.01 V.Ivanchenko remove RogueWave 
//      
// -------------------------------------------------------------------

#include "G4ThirdLevel.hh"


G4ThirdLevel::~G4ThirdLevel(){

  //  this->clearAndDestroy();
  this->clear();
}

G4bool G4ThirdLevel::operator == (const G4ThirdLevel& input) const{

  //  return( this->entries() == input.entries());
  return( this->size() == input.size());

}

G4bool G4ThirdLevel::operator < (const G4ThirdLevel& input) const{

  //  return(this->entries() < input.entries());
  return(this->size() < input.size());

}






