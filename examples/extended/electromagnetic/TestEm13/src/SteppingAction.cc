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
// $Id: SteppingAction.cc,v 1.2 2006-04-14 16:26:43 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "SteppingMessenger.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(PrimaryGeneratorAction* prim,
                               RunAction* RuAct, HistoManager* Hist)
:primary(prim),runAction(RuAct), histoManager(Hist)
{
  stepMessenger = new SteppingMessenger(this);
  fract = 0.1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
  delete stepMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  G4double stepLength = aStep->GetStepLength();

  //count real processes (at least 1 secondary) and sum track length
  //
  G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  if ((*secondary).size() > 0) { 
    runAction->CountProcesses(procName);  
    runAction->SumTrack(stepLength);
  }
  else if (procName == "Transportation") runAction->CountProcesses(procName);
  
  //plot final state (only if continuous energy loss is small enough)
  //
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double E0 = primary->GetParticleGun()->GetParticleEnergy();
  if (edep < fract*E0) {
   
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
    //G4TrackVector* secondary = fpSteppingManager->GetSecondary();
    for (size_t lp=0; lp<(*secondary).size(); lp++) {
      G4double charge = (*secondary)[lp]->GetDefinition()->GetPDGCharge();
      if (charge != 0.) id = 3; else id = 5;
      G4double energy = (*secondary)[lp]->GetKineticEnergy();
      histoManager->FillHisto(id,energy);

      id++;
      G4ThreeVector direction = (*secondary)[lp]->GetMomentumDirection();      
      G4double costeta = direction.x();
      histoManager->FillHisto(id,costeta);         
    }
  }
         
  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


