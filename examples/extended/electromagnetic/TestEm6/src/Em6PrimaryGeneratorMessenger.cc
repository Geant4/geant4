// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em6PrimaryGeneratorMessenger.cc,v 1.1 2000-03-01 07:50:10 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em6PrimaryGeneratorMessenger.hh"

#include "Em6PrimaryGeneratorAction.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6PrimaryGeneratorMessenger::Em6PrimaryGeneratorMessenger(Em6PrimaryGeneratorAction* Em6Gun)
:Em6Action(Em6Gun)
{ 
  DefaultCmd = new G4UIcmdWithoutParameter("/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
  DefaultCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em6PrimaryGeneratorMessenger::~Em6PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em6PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == DefaultCmd )
   { Em6Action->SetDefaultKinematic();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

