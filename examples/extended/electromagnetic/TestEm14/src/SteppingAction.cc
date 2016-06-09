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
// $Id: SteppingAction.cc,v 1.7 2010-10-13 13:42:33 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(PrimaryGeneratorAction* prim,
                               RunAction* RuAct, HistoManager* Hist)
:primary(prim),runAction(RuAct), histoManager(Hist)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();     
  G4bool transmit = (endPoint->GetStepStatus() <= fGeomBoundary);  
  if (transmit) { runAction->CountProcesses(procName); }
  else {                         
    //count real processes and sum track length
    G4double stepLength = aStep->GetStepLength();
    runAction->CountProcesses(procName);  
    runAction->SumTrack(stepLength);
  }
  
  //plot final state (only if continuous energy loss is small enough)
  //
   
  //scattered primary particle
  //
  G4int id = 1;
  if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
    G4double energy = endPoint->GetKineticEnergy();      
    histoManager->FillHisto(id,energy);

    id = 2;
    G4ThreeVector direction = endPoint->GetMomentumDirection();
    G4double costeta = direction.x();
    histoManager->FillHisto(id,costeta);     
  }  
  
  //secondaries
  //
  const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    G4double charge = (*secondary)[lp]->GetDefinition()->GetPDGCharge();
    if (charge != 0.) { id = 3; } else { id = 5; }
    G4double energy = (*secondary)[lp]->GetKineticEnergy();
    histoManager->FillHisto(id,energy);

    ++id;
    G4ThreeVector direction = (*secondary)[lp]->GetMomentumDirection();      
    G4double costeta = direction.x();
    histoManager->FillHisto(id,costeta);
      
    //energy tranferred to charged secondaries
    if (charge != 0.) { runAction->SumeTransf(energy); }         
  }
         
  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


