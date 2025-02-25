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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

    G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()
                            ->GetVolume()->GetName();

    //if a particle enters the detector volume
    if (step->GetPreStepPoint()-> GetStepStatus()==G4StepStatus::fGeomBoundary&&
        (volumeName=="Crystal"||volumeName=="Detector"))
    {

        //coordinates
        G4double x0 = step->GetPreStepPoint()->GetPosition().getX()/CLHEP::mm;
        G4double y0 = step->GetPreStepPoint()->GetPosition().getY()/CLHEP::mm;

        //angles
        G4ThreeVector momentumDirection =
            step->GetPreStepPoint()->GetMomentumDirection();
        G4double angle_x =
            std::atan(momentumDirection.getX()/momentumDirection.getZ());
        G4double angle_y =
            std::atan(momentumDirection.getY()/momentumDirection.getZ());
        if (momentumDirection.getZ() < 0)
        {
            if (momentumDirection.getX() > 0)
            {
                angle_x += CLHEP::pi;
            }
            else
            {
                angle_x -= CLHEP::pi;
            }
            if (momentumDirection.getY() > 0)
            {
                angle_y += CLHEP::pi;
            }
            else
            {
                angle_y -= CLHEP::pi;
            }
        }

        //kinetic energy
        G4double ekin = step->GetPreStepPoint()->GetKineticEnergy()/CLHEP::MeV;

        //particle name, ID and parentID
        G4String particleName =
            step->GetTrack()->GetDefinition()->GetParticleName();
        G4int particleID = step->GetTrack()->GetTrackID();
        G4int parentID = step->GetTrack()->GetParentID();

        //event ID
        G4int eventID =
            G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

        //ntuple index
        G4int iTuple = 0;
        if(volumeName=="Crystal")
        {
            iTuple = 0;
        }
        else if(volumeName=="Detector")
        {
            if(particleName=="gamma")
            {
               iTuple = 2;
            }
            else
            {
               iTuple = 1;
            }
        }

        //saving result to root
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillNtupleIColumn(iTuple,0,eventID);
        analysisManager->FillNtupleSColumn(iTuple,1,volumeName);
        analysisManager->FillNtupleDColumn(iTuple,2,x0);
        analysisManager->FillNtupleDColumn(iTuple,3,y0);
        analysisManager->FillNtupleDColumn(iTuple,4,angle_x);
        analysisManager->FillNtupleDColumn(iTuple,5,angle_y);
        analysisManager->FillNtupleDColumn(iTuple,6,ekin);
        analysisManager->FillNtupleSColumn(iTuple,7,particleName);
        analysisManager->FillNtupleIColumn(iTuple,8,particleID);
        analysisManager->FillNtupleIColumn(iTuple,9,parentID);
        analysisManager->AddNtupleRow(iTuple);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//}
