// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5SteppingMessenger.cc,v 1.1 1999-10-12 12:23:38 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em5SteppingMessenger.hh"

#include "Em5SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5SteppingMessenger::Em5SteppingMessenger(Em5SteppingAction* SA)
:steppingAction (SA)
{
  steppingDir = new G4UIdirectory("/stepping/");
  steppingDir->SetGuidance("stepping control");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5SteppingMessenger::~Em5SteppingMessenger()
{
  delete steppingDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5SteppingMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
