// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: hTestPrimaryGeneratorMessenger.cc,v 1.1 2000-05-21 16:22:19 chauvie Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestPrimaryGeneratorMessenger.hh"

#include "hTestPrimaryGeneratorAction.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorMessenger::hTestPrimaryGeneratorMessenger(hTestPrimaryGeneratorAction* hTestGun)
:hTestAction(hTestGun)
{ 
  DefaultCmd = new G4UIcmdWithoutParameter("/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
  DefaultCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPrimaryGeneratorMessenger::~hTestPrimaryGeneratorMessenger()
{
  delete DefaultCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == DefaultCmd )
   { hTestAction->SetDefaultKinematic();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

