// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07PhysicsListMessenger.cc,v 1.1 1999-01-08 16:35:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07PhysicsListMessenger.hh"

#include "T07PhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07PhysicsListMessenger::T07PhysicsListMessenger(T07PhysicsList * List)
:T07List(List)
{
  ProcessCmd = new G4UIcmdWithoutParameter("/run/setEmProcess",this);
  ProcessCmd->SetGuidance("select electromagnetic processes");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07PhysicsListMessenger::~T07PhysicsListMessenger()
{
  delete ProcessCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void T07PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == ProcessCmd) T07List->SetStatusEmProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

