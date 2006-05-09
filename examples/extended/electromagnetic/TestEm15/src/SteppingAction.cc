//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: SteppingAction.cc,v 1.1 2006-05-09 14:03:04 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det,
                               PrimaryGeneratorAction* prim, RunAction* RuAct, 
			       HistoManager* Hist)
:detector(det), primary(prim), runAction(RuAct), histoManager(Hist)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  
  // if World --> return
  if (prePoint->GetPhysicalVolume() == detector->GetWorld()) return;
  
  // here we enter in the absorber Box
  // kill the event anyway after the first step
  //
  G4RunManager::GetRunManager()->AbortEvent();
  
  //count processes and keep only Multiple Scattering
  //  
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  runAction->CountProcesses(procName);
      
  if (procName != "msc" && procName != "stepMax") return;
  
  //below, only multiple Scattering happens
  //
  G4ThreeVector position  = endPoint->GetPosition();
  G4ThreeVector direction = endPoint->GetMomentumDirection();
  
  G4double truePathLength = aStep->GetStepLength();      
  G4double geomPathLength = position.x() + 0.5*detector->GetBoxSize();
  G4double ratio = geomPathLength/truePathLength;
  runAction->SumPathLength(truePathLength,geomPathLength);
  histoManager->FillHisto(1,truePathLength);
  histoManager->FillHisto(2,geomPathLength);
  histoManager->FillHisto(3,ratio);
   
  G4double yend = position.y(), zend = position.z();
  G4double lateralDisplacement = std::sqrt(yend*yend + zend*zend);
  runAction->SumLateralDisplacement(lateralDisplacement);
  histoManager->FillHisto(4,lateralDisplacement);
  histoManager->FillHisto(5,yend);
  
  G4double psiPlane = std::atan(yend/geomPathLength); 
  runAction->SumPsiPlane(psiPlane);
  histoManager->FillHisto(6,psiPlane);
  
  G4double tetaPlane = std::atan2(direction.y(), direction.x()); 
  runAction->SumTetaPlane(tetaPlane);
  histoManager->FillHisto(7,tetaPlane);
  tetaPlane = std::atan2(direction.z(), direction.x()); 
  runAction->SumTetaPlane(tetaPlane);
  histoManager->FillHisto(7,tetaPlane);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


