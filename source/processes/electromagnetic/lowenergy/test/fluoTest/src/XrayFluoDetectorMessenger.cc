//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorMessenger::XrayFluoDetectorMessenger(XrayFluoDetectorConstruction * Det)
:Detector(Det)
{ 
   detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(Idle);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoDetectorMessenger::~XrayFluoDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


