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
/// \file electromagnetic/TestEm11/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
:G4UserSteppingAction(),fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 G4double edep = step->GetTotalEnergyDeposit();
 if (edep <= 0.) return;
 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();     
 
 //longitudinal profile of deposited energy
 //randomize point of energy deposotion
 //
 G4StepPoint* prePoint  = step->GetPreStepPoint();
 G4StepPoint* postPoint = step->GetPostStepPoint(); 
 G4ThreeVector P1 = prePoint ->GetPosition();
 G4ThreeVector P2 = postPoint->GetPosition();
 G4ThreeVector point = P1 + G4UniformRand()*(P2 - P1);
 if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0.) point = P2;
 G4double x = point.x();
 G4double xshifted = x + 0.5*fDetector->GetAbsorSizeX();
 analysisManager->FillH1(1, xshifted, edep);

 //"normalized" histogram
 // 
 Run* run
   = static_cast<Run*>(
       G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 G4int iabs = prePoint->GetTouchableHandle()->GetCopyNumber(1);
 G4double csdaRange  = run->GetCsdaRange(iabs);
 if (csdaRange > 0.) { 
   G4double density = fDetector->GetAbsorMaterial(iabs)->GetDensity();
   G4double xfront  = fDetector->GetXfront(iabs);
   G4double xfrontNorm = run->GetXfrontNorm(iabs);
   G4double xnorm = xfrontNorm + (x - xfront)/csdaRange;
   analysisManager->FillH1(8, xnorm, edep/(csdaRange*density));
 }
 
 //total energy deposit in absorber
 //
 fEventAction->AddEdep(iabs, edep);
 
 //step size of primary particle or charged secondaries
 //
 G4double steplen = step->GetStepLength();
 const G4Track* track = step->GetTrack();
 if      (track->GetTrackID() == 1) analysisManager->FillH1(4, steplen);
 else if (track->GetDefinition()->GetPDGCharge() != 0.)
                                    analysisManager->FillH1(7, steplen);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

