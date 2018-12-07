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
/// \file TrackingAction.hh
/// \brief Implementation of the TrackingAction class.

#include "TrackingAction.hh"
#include "UserTrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

TrackingAction::TrackingAction(): G4UserTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    // All primaries have a flag/label equal to 1
    if ( aTrack->GetParentID() == 0 ) {
        UserTrackInformation* parentInformation = new UserTrackInformation();
        parentInformation->SetSplitTrackID(1);
        G4Track* newTrack = (G4Track*)aTrack;
        newTrack->SetUserInformation(parentInformation);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
    UserTrackInformation* parentInformation = 
                                 (UserTrackInformation*)(aTrack->GetUserInformation());
    // Pass the flag to all seconaries
    G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
    if (secondaries)
    {
        size_t nbSec = secondaries->size();
        for (size_t iSec=0; iSec < nbSec; iSec++) {
            if ((*secondaries)[iSec]->GetUserInformation() == 0 ) {
                // Retrieve the flag 
                UserTrackInformation* secondaryInformation = new UserTrackInformation();
                (*secondaries)[iSec]->SetUserInformation(secondaryInformation);
                
                G4int splitTrackID = parentInformation->GetSplitTrackID();
                // Pass the flag
                secondaryInformation->SetSplitTrackID(splitTrackID);
            }
        }
    }
}

void TrackingAction::Initialize() {
        return;
}
