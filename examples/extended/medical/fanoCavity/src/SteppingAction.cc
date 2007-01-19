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
// $Id: SteppingAction.cc,v 1.1 2007-01-19 17:20:27 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct,
                               HistoManager* histo)
:detector(det), runAction(RuAct), histoManager(histo),
 wall(0), cavity(0)
{ 
  first = true;
  trackCharged = false;
  trackSegm = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
 //get detector pointers
 if (first) {
   wall   = detector->GetWall();
   cavity = detector->GetCavity();
   first  = false;
 }
 
 //get volume
 //
 G4StepPoint* point1 = step->GetPreStepPoint();
 G4VPhysicalVolume* volume = point1->GetTouchableHandle()->GetVolume();
  
 // count processes
 //
 G4StepPoint* point2 = step->GetPostStepPoint(); 
 const G4VProcess* process = point2->GetProcessDefinedStep();
 if (process) runAction->CountProcesses(process->GetProcessName());
 
 //energy deposit in cavity
 //
 if (volume == cavity) { 
   G4double edep = step->GetTotalEnergyDeposit();
   if (edep > 0.) runAction->AddEdepCavity(edep);     
 }

 //keep only charged particles
 //
 if (!trackCharged) return;
 
 //step size of charged particles
 //
 G4int id;
 G4double steplen = step->GetStepLength();
 if (volume == wall) {runAction->StepInWall  (steplen); id = 8;}
 else                {runAction->StepInCavity(steplen); id = 9;}
 histoManager->FillHisto(id,steplen);
 
 //keep only charged particles within cavity
 //
 if (volume == wall) return;
 
 G4double ekin1 = point1->GetKineticEnergy();
 G4double ekin2 = point2->GetKineticEnergy();
 
 //enter in cavity
 //
 if (point1->GetStepStatus() == fGeomBoundary) {
   runAction->FlowInCavity(0,ekin1);
   histoManager->FillHisto(5,ekin1);    
   G4ThreeVector direction = point1->GetMomentumDirection();
   histoManager->FillHisto(6,std::acos(direction.z()));  
   G4ThreeVector vertex = step->GetTrack()->GetVertexPosition();
   histoManager->FillHisto(4,vertex.z());   
   trackSegm = 0.;      
 }
  
 //within cavity
 //
 if (step->GetTrack()->GetCurrentStepNumber() == 1) trackSegm = 0.;
 trackSegm += steplen;
 if (ekin2 <= 0.) {
   runAction->AddTrakCavity(trackSegm);
   histoManager->FillHisto(7,trackSegm);      
 } 
 
 //exit cavity
 //
 if (point2->GetStepStatus() == fGeomBoundary) {
   runAction->FlowInCavity(1,ekin2);
   runAction->AddTrakCavity(trackSegm);
   histoManager->FillHisto(7,trackSegm);      
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


