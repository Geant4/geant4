// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4SteppingAction.cc,v 1.1 1999-10-12 11:26:59 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em4SteppingAction.hh"
#include "Em4EventAction.hh"
#include "G4SteppingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4SteppingAction::Em4SteppingAction(Em4EventAction* EvAct)
:eventAction(EvAct),myVerbose(NULL)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em4SteppingAction::~Em4SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em4SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 G4double EdepStep = aStep->GetTotalEnergyDeposit();
 if (EdepStep > 0.) eventAction->addEdep(EdepStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


