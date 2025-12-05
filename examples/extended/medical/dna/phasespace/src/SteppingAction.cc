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
/// \brief Implementation of the SteppingAction class

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "SteppingAction.hh"
#include "SteppingMessenger.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4Alpha.hh"
#include "G4AnalysisManager.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det) : G4UserSteppingAction(), 
fpDetectorConstruction(det)
{
  fpSteppingMessenger = new SteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() 
{
  delete fpSteppingMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Protection
  if (!step->GetPostStepPoint()) return;
  //
  
  // PreStep point
  
  G4StepPoint* preStep = step->GetPreStepPoint();

  const G4LogicalVolume* preVolume =
    preStep->GetPhysicalVolume()->GetLogicalVolume();

  // PostStep point
  
  G4StepPoint* postStep = step->GetPostStepPoint();
 
  // Protection
  if (postStep->GetStepStatus() == fWorldBoundary ) return;
  //

  const G4LogicalVolume* postVolume = 
    postStep->GetPhysicalVolume()->GetLogicalVolume();

  // METHOD 1: kill the particle and sec. if inside the scorer
 
   if ( postVolume == fpDetectorConstruction->GetLogicalScorer() &&
        preVolume  == fpDetectorConstruction->GetLogicalScorer() )
  {    
    if (fKill) step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
  }

  // Information is scored from preStep = WORLD to postStep = SCORER 

  if ( preStep->GetStepStatus() == fGeomBoundary )
  {    
    // G4cout << "---> One particle has crossed the boundary..." << G4endl;

    G4double PDGCode = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  
    G4double x = preStep->GetPosition().x() / micrometer;
    G4double y = preStep->GetPosition().y() / micrometer;
    G4double z = preStep->GetPosition().z() / micrometer;

    G4ThreeVector direction = preStep->GetMomentumDirection();
    G4double xMom = direction.x();
    G4double yMom = direction.y();
    G4double zMom = direction.z();
  
    G4double energy = preStep->GetKineticEnergy() / keV;

    //
  
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    analysisManager->FillNtupleDColumn(0, PDGCode);
    analysisManager->FillNtupleDColumn(1, x);
    analysisManager->FillNtupleDColumn(2, y);
    analysisManager->FillNtupleDColumn(3, z);
    analysisManager->FillNtupleDColumn(4, xMom);
    analysisManager->FillNtupleDColumn(5, yMom);
    analysisManager->FillNtupleDColumn(6, zMom);
    analysisManager->FillNtupleDColumn(7, energy);

    //
    analysisManager->AddNtupleRow();

    // METHOD 2: kill the particle and sec. if entering
    // if (fKill) step->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);

  } // End of scoring

}
