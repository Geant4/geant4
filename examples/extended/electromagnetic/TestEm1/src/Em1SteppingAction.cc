// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1SteppingAction.cc,v 1.2 1999-12-15 14:48:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1SteppingAction.hh"
#include "Em1RunAction.hh"
#include "Em1EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"

#include "CLHEP/Hist/HBookFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1SteppingAction::Em1SteppingAction(Em1RunAction* RuAct, Em1EventAction* EvAct)
:runAction(RuAct),eventAction(EvAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1SteppingAction::~Em1SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 
 G4double EdepStep = aStep->GetTotalEnergyDeposit();
 if (EdepStep > 0.) eventAction->addEdep(EdepStep);
 
 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
 if (process) runAction->CountProcesses(process->GetProcessName());
 
 G4double charge  = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
 G4double steplen = aStep->GetStepLength();
 if (charge != 0.) runAction->GetHisto(2)->accumulate(steplen);                   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


