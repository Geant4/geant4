
#include "Tst10DetectorMessenger.hh"

#include "Tst10DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

Tst10DetectorMessenger::Tst10DetectorMessenger(Tst10DetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selDetCmd = new G4UIcmdWithAString("/mydet/SelectDetector",this);
  selDetCmd->SetGuidance("Select Detector Setup.");
  selDetCmd->SetGuidance("  Choice : Detector type ");
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("Sphere");
  selDetCmd->SetCandidates("Sphere Box Cone Tube Hype Torus Para Trd");
  selDetCmd->AvailableForStates(PreInit,Idle);

  myDetector->SelectDetector(defParam="Sphere");
}

void Tst10DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
		myDetector->SwitchDetector();
  }
  return;
}

