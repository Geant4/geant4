//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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

  G4int idl=0, idn=0;
  G4int level = theTouchable->GetHistoryDepth();  
  G4int idx = theTouchable->GetReplicaNumber( 0 ) - 1;
  if ( level > 0 ) {
    idl = theTouchable->GetReplicaNumber( 1 ) - 1;
    if ( level > 1 ) {
      idn = theTouchable->GetReplicaNumber( 2 );
    }
  }

  unsigned int idunit = idn*4096 + idl*64 + idx;

  return idunit;
}
