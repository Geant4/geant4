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
//
// $Id: SteppingAction.cc,v 1.1 2004/04/28 11:12:40 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

#ifdef USE_AIDA
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct, HistoManager* Hist)
  :runAction(RuAct), histoManager(Hist)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  G4double stepLength = aStep->GetStepLength();

  //count processes and sum track length
  //
  runAction->CountProcesses(procName);
  if (procName != "Transportation") runAction->SumTrack(stepLength);
  
#ifdef USE_AIDA  
  //plot final state
  //
  //scattered primary particle
  //
  G4int id = 1;
  if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
    if (histoManager->GetHisto(id)) {
      G4double energy = endPoint->GetKineticEnergy();      
      G4double unit = histoManager->GetHistoUnit(id); 
      histoManager->GetHisto(id)->fill(energy/unit);
    }
    id = 2;
    if (histoManager->GetHisto(id)) {
      G4ThreeVector direction = endPoint->GetMomentumDirection();
      G4double costeta = direction.x();
      histoManager->GetHisto(id)->fill(costeta); 
    }        
  }  
  
  //secondaries
  //
  G4TrackVector* secondary = fpSteppingManager->GetSecondary();
  for (size_t lp=0; lp<(*secondary).size(); lp++) {
    G4double charge = (*secondary)[lp]->GetDefinition()->GetPDGCharge();
    if (charge != 0.) id = 3; else id = 5;
    if (histoManager->GetHisto(id)) {
      G4double energy = (*secondary)[lp]->GetKineticEnergy();
      G4double unit = histoManager->GetHistoUnit(id);              
      histoManager->GetHisto(id)->fill(energy/unit);
    }
    id++;
    if (histoManager->GetHisto(id)) {
      G4ThreeVector direction = (*secondary)[lp]->GetMomentumDirection();      
      G4double costeta = direction.x();
      histoManager->GetHisto(id)->fill(costeta);    
    }                       
  }
#endif
         
  // kill event after first interaction
  //
  G4RunManager::GetRunManager()->AbortEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


