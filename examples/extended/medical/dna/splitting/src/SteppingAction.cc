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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class.

#include "SteppingAction.hh"
#include "UserTrackInformation.hh"

#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction(),fpEventAction(0)
{
    fRegion = G4RegionStore::GetInstance()->FindOrCreateRegion("Target");
    fpEventAction = (EventAction*) G4EventManager::GetEventManager()->
                  GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
    if(theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()!=
        "Transportation")
    {
        // Get if hit was an ionization and get energy deposit
        //
        G4double edepStep = theStep->GetTotalEnergyDeposit()/eV;
        G4Region* currentRegion = 
            theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->
            GetLogicalVolume()->GetRegion();

        if (edepStep > 0. && (fRegion == currentRegion ) ) 
        {
            const G4String processName = theStep->GetPostStepPoint()->
                                         GetProcessDefinedStep()->GetProcessName();

            if ( "Wrappede-_G4DNAIonisation" == processName || 
                 "e-_G4DNAIonisation" == processName ) {  
                // Get the flag/label
                UserTrackInformation* trackInformation = 
                      (UserTrackInformation*)(theStep->GetTrack()->GetUserInformation());
                G4int idx = trackInformation->GetSplitTrackID();
               
                // score hits for each flag
                fpEventAction->AddIonizationEvent(idx, 1);
                fpEventAction->AddEdepEvent(
                               edepStep*theStep->GetPreStepPoint()->GetWeight());
            }
        }
    }
}

