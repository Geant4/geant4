//  XrayTelPrimaryGeneratorMessenger.cc
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayTelPrimaryGeneratorMessenger.hh"

#include "XrayTelPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorMessenger::XrayTelPrimaryGeneratorMessenger(XrayTelPrimaryGeneratorAction* XrayTelGun)
:XrayTelAction(XrayTelGun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : isotropic, beam, aperture, point");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("isotropic");
  RndmCmd->SetCandidates("isotropic beam aperture point");
  RndmCmd->AvailableForStates(PreInit,Idle);

  ErndmCmd = new G4UIcmdWithAString("/gun/erandom",this);
  ErndmCmd->SetGuidance("Select particle energy at random");
  ErndmCmd->SetGuidance("  Choice : rndm mono");
  ErndmCmd->SetParameterName("choice",true);
  ErndmCmd->SetDefaultValue("mono");
  ErndmCmd->SetCandidates(" rndm mono");
  ErndmCmd->AvailableForStates(PreInit,Idle);

  SetRmin  = new G4UIcmdWithADoubleAndUnit("/gun/SetRmin",this);
  SetRmin->SetGuidance("Set the minimum radius for random position generation");
  SetRmin->SetParameterName("Rmin",true);
  SetRmin->SetDefaultUnit("cm"); 
  SetRmin->SetDefaultValue(0.0*cm);

  SetRmax = new G4UIcmdWithADoubleAndUnit("/gun/SetRmax",this);
  SetRmax->SetGuidance("Set the maximum radius for random position generation");
  SetRmax->SetParameterName("Rmax",true);
  SetRmax->SetDefaultUnit("cm"); 
  SetRmax->SetDefaultValue(100.0*cm);

  SetTmax = new G4UIcmdWithADoubleAndUnit("/gun/SetTmax",this);
  SetTmax->SetGuidance("Set the maximum aperture angle for direction generation");
  SetTmax->SetParameterName("Tmax",true);
  SetTmax->SetDefaultUnit("degree"); 
  SetTmax->SetDefaultValue(1*degree);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelPrimaryGeneratorMessenger::~XrayTelPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete ErndmCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == RndmCmd ){ 
    XrayTelAction->SetRndmFlag(newValue);
  } else if( command == ErndmCmd ){ 
    XrayTelAction->SetErndmFlag(newValue);  
  } else if  ( command == SetRmin ){
    XrayTelAction->SetRmin( SetRmin->GetNewDoubleValue ( newValue )  );
  } else if  ( command == SetRmax ){
    XrayTelAction->SetRmax( SetRmax->GetNewDoubleValue ( newValue )  );
  } else if ( command == SetTmax ){
    XrayTelAction->SetTmax( SetTmax->GetNewDoubleValue ( newValue )  );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






