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
///////////////////////////////////////////////////////////////////////////////
// File: CCalEcalOrganization.cc
// Description: Defines numbering schema for the Electromagnetic Calorimeter
///////////////////////////////////////////////////////////////////////////////
#include "CCalEcalOrganization.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"

CCalEcalOrganization::~CCalEcalOrganization(){
  G4cout << " Deleting CCalEcalOrganization" << G4endl;
}


unsigned int CCalEcalOrganization::GetUnitID(const G4Step* aStep) const {

  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );

  int idl=0, idn=0;
  int level = theTouchable->GetHistoryDepth();  
  int idx = theTouchable->GetReplicaNumber( 0 ) - 1;
  if ( level > 0 ) {
    idl = theTouchable->GetReplicaNumber( 1 ) - 1;
    if ( level > 1 ) {
      idn = theTouchable->GetReplicaNumber( 2 );
    }
  }

  unsigned int idunit = idn*4096 + idl*64 + idx;

  return idunit;
}
