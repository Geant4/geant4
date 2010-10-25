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
// $Id: SteppingAction.cc,v 1.3 2010-10-25 13:31:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "HistoManager.hh"

#include "G4SteppingManager.hh"
#include "G4Gamma.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* RuAct,
                               EventAction* EvAct,TrackingAction* TrAct,
			       HistoManager* histo)
:detector(det), runAction(RuAct), eventAction(EvAct), trackAction(TrAct),
 histoManager(histo), wall(0), cavity(0)
{ 
  first = true;
  trackSegm = 0.;
  directionIn = G4ThreeVector(0.,0.,0.);
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
   if (edep > 0.) trackAction->AddEdepCavity(edep);     
 }

 //keep only charged particles
 //
 if (step->GetTrack()->GetDefinition() == G4Gamma::Gamma()) return;
 
 //step size of charged particles
 //
 G4int id;
 G4double steplen = step->GetStepLength();
 if (volume == wall) {runAction->StepInWall  (steplen); id = 9;}
 else                {runAction->StepInCavity(steplen); id = 10;}
 histoManager->FillHisto(id,steplen);
 
 //last step before hitting the cavity
 //
 if ((volume == wall) && (point2->GetStepStatus() == fGeomBoundary)) {
   directionIn = point1->GetMomentumDirection();
 } 
 
 //keep only charged particles within cavity
 //
 if (volume == wall) return;
 
 G4double ekin1 = point1->GetKineticEnergy();
 G4double ekin2 = point2->GetKineticEnergy();
 
 //first step in cavity
 //
 if (point1->GetStepStatus() == fGeomBoundary) {
   trackSegm = 0.;
   G4ThreeVector vertex = step->GetTrack()->GetVertexPosition();
   histoManager->FillHisto(4,vertex.z());          
   runAction->FlowInCavity(0,ekin1);
   histoManager->FillHisto(5,ekin1);
   if (steplen>0.) {    
     G4ThreeVector directionOut = 
              (point2->GetPosition() - point1->GetPosition()).unit();
     G4ThreeVector normal = point1->GetTouchableHandle()->GetSolid()
                            ->SurfaceNormal(point1->GetPosition());
     histoManager->FillHisto(6,std::acos(-directionIn*normal));
     histoManager->FillHisto(7,std::acos(-directionOut*normal));
   }		   
 }
  
 //within cavity
 //
 if (step->GetTrack()->GetCurrentStepNumber() == 1) trackSegm = 0.;
 trackSegm += steplen;
 if (ekin2 <= 0.) {
   runAction->AddTrakCavity(trackSegm);
   histoManager->FillHisto(8,trackSegm);      
 } 
 
 //exit cavity
 //
 if (point2->GetStepStatus() == fGeomBoundary) {
   runAction->FlowInCavity(1,ekin2);
   runAction->AddTrakCavity(trackSegm);
   histoManager->FillHisto(8,trackSegm);      
 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


