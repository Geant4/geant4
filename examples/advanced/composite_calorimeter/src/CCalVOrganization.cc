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
// File: CCalVOrganization.cc
// Description: Base class for definition of sensitive unit numbering schema
///////////////////////////////////////////////////////////////////////////////
#include "CCalVOrganization.hh"
#include "G4TouchableHistory.hh"
#include "globals.hh"

G4int CCalVOrganization::Levels(const G4Step* aStep) const {

  //Find number of levels
  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  return ((theTouchable->GetHistoryDepth()) + 1);
}


void CCalVOrganization::DetectorLevel(const G4Step* aStep, G4int& level,
                                     G4int* copyno, G4String* name) const {

  //Get name and copy numbers
  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  for ( G4int ii = 0; ii < level; ++ii ) {
    name[ level - ii - 1 ]   = theTouchable->GetVolume(ii)->GetName();
    copyno[ level - ii - 1 ] = theTouchable->GetReplicaNumber(ii);
  }

}
