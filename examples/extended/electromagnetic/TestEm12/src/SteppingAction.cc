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
/// \file electromagnetic/TestEm12/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* event)
:G4UserSteppingAction(),
 fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 G4double edep = step->GetTotalEnergyDeposit();
 if (edep <= 0.) return;
 
 //total energy deposit in absorber
 //
 fEventAction->AddEdep(edep);     
 
 //longitudinal profile of deposited energy
 //        
 G4ThreeVector prePoint  = step->GetPreStepPoint() ->GetPosition();
 G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
 G4ThreeVector point = prePoint + G4UniformRand()*(postPoint - prePoint);
 if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0) point = postPoint;
 G4double r = point.mag();
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->FillH1(1, r, edep);
 
 Run* run
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    
 G4double r0 = run->GetCsdaRange();
 if (r0 > 0.) analysisManager->FillH1(8, r/r0, edep);
 
 //step size of primary particle or charged secondaries
 //
 G4double steplen = step->GetStepLength();
 const G4Track* track = step->GetTrack();
 if      (track->GetTrackID() == 1) analysisManager->FillH1(4, steplen);
 else if (track->GetDefinition()->GetPDGCharge() != 0.)
                                    analysisManager->FillH1(7, steplen); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


