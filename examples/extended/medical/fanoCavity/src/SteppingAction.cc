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
/// \file medical/fanoCavity/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "TrackingAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "G4RunManager.hh"

#include "G4SteppingManager.hh"
#include "G4Gamma.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,TrackingAction* TrAct)
                               :fDetector(det), fTrackAction(TrAct),
                                fWall(0), fCavity(0)
{ 
  first = true;
  fTrackSegm = 0.;
  fDirectionIn = G4ThreeVector(0.,0.,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 //get fDetector pointers
 if (first) {
   fWall   = fDetector->GetWall();
   fCavity = fDetector->GetCavity();
   first  = false;
 }

 //histograms 
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
 //get volume
 //
 G4StepPoint* point1 = step->GetPreStepPoint();
 G4VPhysicalVolume* volume = point1->GetTouchableHandle()->GetVolume();
  
 Run* run = static_cast<Run*>(
            G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 // count processes
 //
 G4StepPoint* point2 = step->GetPostStepPoint(); 
 const G4VProcess* process = point2->GetProcessDefinedStep();
 if (process) run->CountProcesses(process->GetProcessName());
 
 //energy deposit in cavity
 //
 if (volume == fCavity) { 
   G4double edep = step->GetTotalEnergyDeposit();
   if (edep > 0.) fTrackAction->AddEdepCavity(edep);     
 }

 //keep only charged particles
 //
 if (step->GetTrack()->GetDefinition() == G4Gamma::Gamma()) return;
 
 //step size of charged particles
 //
 G4int id;
 G4double steplen = step->GetStepLength();
 if (volume == fWall) {run->StepInWall  (steplen); id = 9;}
 else                {run->StepInCavity(steplen); id = 10;}
 analysisManager->FillH1(id,steplen);
 
 //last step before hitting the cavity
 //
 if ((volume == fWall) && (point2->GetStepStatus() == fGeomBoundary)) {
   fDirectionIn = point1->GetMomentumDirection();
 } 
 
 //keep only charged particles within cavity
 //
 if (volume == fWall) return;
 
 G4double ekin1 = point1->GetKineticEnergy();
 G4double ekin2 = point2->GetKineticEnergy();
 
 //first step in cavity
 //
 if (point1->GetStepStatus() == fGeomBoundary) {
   fTrackSegm = 0.;
   G4ThreeVector vertex = step->GetTrack()->GetVertexPosition();
   analysisManager->FillH1(4,vertex.z());          
   run->FlowInCavity(0,ekin1);
   analysisManager->FillH1(5,ekin1);
   if (steplen>0.) {    
     G4ThreeVector directionOut = 
              (point2->GetPosition() - point1->GetPosition()).unit();
     G4ThreeVector normal = point1->GetTouchableHandle()->GetSolid()
                            ->SurfaceNormal(point1->GetPosition());
     analysisManager->FillH1(6,std::acos(-fDirectionIn*normal));
     analysisManager->FillH1(7,std::acos(-directionOut*normal));
   }                   
 }
  
 //within cavity
 //
 if (step->GetTrack()->GetCurrentStepNumber() == 1) fTrackSegm = 0.;
 fTrackSegm += steplen;
 if (ekin2 <= 0.) {
   run->AddTrakCavity(fTrackSegm);
   analysisManager->FillH1(8,fTrackSegm);      
 } 
 
 //exit cavity
 //
 if (point2->GetStepStatus() == fGeomBoundary) {
   run->FlowInCavity(1,ekin2);
   run->AddTrakCavity(fTrackSegm);
   analysisManager->FillH1(8,fTrackSegm);      
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

