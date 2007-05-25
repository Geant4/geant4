#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                             PrimaryGeneratorAction* Gun)
:Action(Gun)
{
  gunDir = new G4UIdirectory("/testem/gun/");
  gunDir->SetGuidance("gun control");
   
  DefaultCmd = new G4UIcmdWithoutParameter("/testem/gun/setDefault",this);
  DefaultCmd->SetGuidance("set/reset kinematic defined in PrimaryGenerator");
  DefaultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete DefaultCmd;
  delete gunDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
					    G4String)
{ 
  if( command == DefaultCmd )
   { Action->SetDefaultKinematic();}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

