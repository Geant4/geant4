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
//
/// \file field/field04/src/F04TrackingAction.cc
/// \brief Implementation of the F04TrackingAction class
//

#include "globals.hh"
#include "G4RunManager.hh"

#include "F04UserTrackInformation.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4TrackingManager.hh"

#include "F04TrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  F04UserTrackInformation* trackInformation = new F04UserTrackInformation();

  fpTrackingManager->SetUserTrackInformation(trackInformation);

  if (aTrack->GetMomentumDirection().z()>0.0) {
     trackInformation->SetTrackStatusFlag(right);
  } else {
     trackInformation->SetTrackStatusFlag(left);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04TrackingAction::PostUserTrackingAction(const G4Track* aTrack){

  F04UserTrackInformation*
    trackInformation = (F04UserTrackInformation*)aTrack->GetUserInformation();
 
  if ( aTrack->GetDefinition()==G4MuonPlus::MuonPlusDefinition() ||
       aTrack->GetDefinition()==G4PionPlus::PionPlusDefinition() ) {
    if (trackInformation->GetTrackStatusFlag() == reverse) {
//       G4RunManager::GetRunManager()->rndmSaveThisEvent();
    }
  }

}
