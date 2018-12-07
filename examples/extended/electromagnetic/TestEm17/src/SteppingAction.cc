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
/// \file electromagnetic/TestEm17/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(RunAction* RuAct, HistoManager* Hist)
  :G4UserSteppingAction(),fRunAction(RuAct), fHistoManager(Hist)
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
  fRunAction->CountProcesses(procName);
  
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
  if (procName == "muIoni")          id = 1; 
  else if (procName == "muPairProd") id = 2;
  else if (procName == "muBrems")    id = 3;
  //  else if (procName == "muNucl")     id = 4;    
  else if (procName == "muonNuclear")id = 4;    
  else if (procName == "hIoni")      id = 5; 
  else if (procName == "hPairProd")  id = 6;
  else if (procName == "hBrems")     id = 7;
  fHistoManager->FillHisto(id,lgepsE);                       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


