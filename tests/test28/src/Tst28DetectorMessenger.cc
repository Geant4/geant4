
#include "Tst28DetectorMessenger.hh"

#include "Tst28DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

Tst28DetectorMessenger::Tst28DetectorMessenger(Tst28DetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selMatCmd = new G4UIcmdWithAString("/mydet/SelectMaterial",this);
  selMatCmd->SetGuidance("Select Material of the SimpleBox.");
  selMatCmd->SetGuidance("  Choice : H, Si, Cu, H2O(Water), U (default)");
  selMatCmd->SetParameterName("choice",true);
  selMatCmd->SetDefaultValue("U");
  selMatCmd->SetCandidates("H Si Cu Pb U H2O");
  selMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  myDetector->SelectMaterial(defParam="U");
}

void Tst28DetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selMatCmd )
  {
    myDetector->SelectMaterial(newValues);
  }
  return;
}

