//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "fluoTestPrimaryGeneratorMessenger.hh"

#include "fluoTestPrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestPrimaryGeneratorMessenger::fluoTestPrimaryGeneratorMessenger(fluoTestPrimaryGeneratorAction* Gun)
:Action(Gun)
{ 
  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : off(default), on");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);
  
  RndmCmmd = new G4UIcmdWithAString("/gun/randomEnergy",this);
  RndmCmmd->SetGuidance("Shoot randomly the energy of the incident particle.");
  RndmCmmd->SetGuidance("  Choice : off(default), on");
  RndmCmmd->SetParameterName("choice",true);
  RndmCmmd->SetDefaultValue("off");
  RndmCmmd->SetCandidates("on off");
  RndmCmmd->AvailableForStates(PreInit,Idle);

  SigmAngleCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaAngle",this);
  SigmAngleCmd->SetGuidance("Select standard deviation of the direction of the     beam.");
  SigmAngleCmd->SetParameterName("choice",false);
  SigmAngleCmd->SetRange("Size>=0.,Size<(2*pi)");
  SigmAngleCmd->SetUnitCategory("Length");
  SigmAngleCmd->AvailableForStates(Idle);
 
  SigmaMomentumCmd = new G4UIcmdWithADoubleAndUnit("/gun/sigmaMomentum",this);
  SigmaMomentumCmd->SetGuidance("Select standard deviation of the momentum of the particles.");
  SigmaMomentumCmd->SetParameterName("choice",false);
  SigmaMomentumCmd->SetRange("Size>=0.");
  SigmaMomentumCmd->SetUnitCategory("Length");
  SigmaMomentumCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

fluoTestPrimaryGeneratorMessenger::~fluoTestPrimaryGeneratorMessenger()
{
  delete RndmCmd;
  delete RndmCmmd;
  delete SigmaMomentumCmd;
  delete SigmAngleCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void fluoTestPrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue )
{ 
  if( command == RndmCmd )
   { Action->SetRndmFlag(newValue);}
 
  if( command == RndmCmmd )
    { Action->SetRandomizePrimary(newValue);}
  if(command == SigmAngleCmd)
    { Action->SetSigmaAngle(SigmAngleCmd->GetNewDoubleValue(newValue));}
  if(command == SigmaMomentumCmd) 
    { Action->SetSigmaMomentum(SigmaMomentumCmd->
     GetNewDoubleValue(newValue));}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...



