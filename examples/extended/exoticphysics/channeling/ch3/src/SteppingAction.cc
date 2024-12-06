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
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

    G4String volumeName = step->GetPreStepPoint()->GetTouchableHandle()
                            ->GetVolume()->GetName();

    //if a particle enters the detector volume
    if (step->GetPreStepPoint()-> GetStepStatus()==G4StepStatus::fGeomBoundary&&
        volumeName=="Detector")
    {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

        //e- spectrum
        if(step->GetTrack()->GetDefinition()->GetParticleName()=="e-")
        {
            //coordinate in the horizontal plane
            G4double eElectron = step->GetPreStepPoint()->GetTotalEnergy()/CLHEP::GeV;

            //filling histogram
            analysisManager->FillH1(0, eElectron);
        }

        //e+ spectrum
        if(step->GetTrack()->GetDefinition()->GetParticleName()=="e+")
        {
            //coordinate in the horizontal plane
            G4double ePositron = step->GetPreStepPoint()->GetTotalEnergy()/CLHEP::GeV;

            //filling histogram
            analysisManager->FillH1(1, ePositron);
        }

        //gamma spectrum
        if(step->GetTrack()->GetDefinition()->GetParticleName()=="gamma")
        {
            //coordinate in the horizontal plane
            G4double eGamma = step->GetPreStepPoint()->GetTotalEnergy()/CLHEP::GeV;

            //filling histogram
            analysisManager->FillH1(2, eGamma);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//}
