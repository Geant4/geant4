// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3PrimaryGeneratorMessenger.cc,v 1.1 1999-10-11 16:55:54 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3PrimaryGeneratorMessenger.hh"

#include "Em3PrimaryGeneratorAction.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3PrimaryGeneratorMessenger::Em3PrimaryGeneratorMessenger(Em3PrimaryGeneratorAction* Em3Gun)
:Em3Action(Em3Gun)
{ 
  DefaultCmd = new G4UIcmdWithoutParameter("/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset the kinematic defined in PrimaryGenerator");
  DefaultCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3PrimaryGeneratorMessenger::~Em3PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == DefaultCmd )
   { Em3Action->SetDefaultKinematic();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

