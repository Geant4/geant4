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
// $Id: SteppingAction.cc,v 1.4 2004/12/03 09:38:31 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct, HistoManager* Hist)
  :runAction(RuAct), histoManager(Hist)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{  
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();

  //count processes 
  //    
  runAction->CountProcesses(procName);
  
  //plot energy transfered
  //
  G4StepPoint* startPoint = aStep->GetPreStepPoint();
  G4double E0 = startPoint->GetKineticEnergy();
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double E1 = E0 - edep;
  G4double E2 = endPoint->GetKineticEnergy();
  G4double etrans = E1 - E2;
  G4double lgepsE = 0.;
  if (etrans > 0.) lgepsE = std::log10(etrans/E1);

  G4int id = 0;
  if (procName == "muIoni")     id = 1; 
  if (procName == "muPairProd") id = 2;
  if (procName == "muBrems")    id = 3;
  if (procName == "muNucl")     id = 4;    
  histoManager->FillHisto(id,lgepsE);		       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


