//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestPhysicsListMessenger.hh"
#include "FluoTestPhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsListMessenger::FluoTestPhysicsListMessenger(FluoTestPhysicsList * List)
:FluoTestList(List)
{

  EnDir = new G4UIdirectory("/energy/");
  EnDir->SetGuidance("energy commands");

  cutGLowLimCmd = new G4UIcmdWithADoubleAndUnit("/energy/lowlimG",this);
  cutGLowLimCmd->SetGuidance("Set ENERGY low limit for Gamma.");
  cutGLowLimCmd->SetParameterName("energy",true);
  cutGLowLimCmd->SetDefaultValue(1e-3);
  cutGLowLimCmd->SetDefaultUnit("MeV");
  cutGLowLimCmd->AvailableForStates(Idle);

  cutELowLimCmd = new G4UIcmdWithADoubleAndUnit("/energy/lowlimE",this);
  cutELowLimCmd->SetGuidance("Set ENERGY low limit for e-.");
  cutELowLimCmd->SetParameterName("energy",true);
  cutELowLimCmd->SetDefaultValue(1e-3);
  cutELowLimCmd->SetDefaultUnit("MeV");
  cutELowLimCmd->AvailableForStates(Idle);

  cutGELowLimCmd = new G4UIcmdWithADoubleAndUnit("/energy/lowlimGE",this);
  cutGELowLimCmd->SetGuidance("Set ENERGY low limit for e- and Gamma.");
  cutGELowLimCmd->SetParameterName("energy",true);
  cutGELowLimCmd->SetDefaultValue(1e-3);
  cutGELowLimCmd->SetDefaultUnit("MeV");
  cutGELowLimCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsListMessenger::~FluoTestPhysicsListMessenger()
{

  delete cutGLowLimCmd;
  delete cutELowLimCmd;
  delete cutGELowLimCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void FluoTestPhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if(command == cutGLowLimCmd)
    { FluoTestList->SetGammaLowLimit(cutGLowLimCmd->GetNewDoubleValue(newValue));}
  if(command == cutELowLimCmd)
    { FluoTestList->SetElectronLowLimit(cutELowLimCmd->GetNewDoubleValue(newValue));}

  if(command == cutGELowLimCmd)
    { FluoTestList->SetGELowLimit(cutGELowLimCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

