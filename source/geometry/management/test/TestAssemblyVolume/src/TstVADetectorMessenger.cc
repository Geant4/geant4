
#include "TstVADetectorMessenger.hh"

#include "TstVADetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

TstVADetectorMessenger::TstVADetectorMessenger(TstVADetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/geom/");
  mydetDir->SetGuidance("Geometry setup commands.");

  selDetCmd = new G4UIcmdWithAString("/geom/select",this);
  selDetCmd->SetGuidance("Select the way detector geometry is built.");
  selDetCmd->SetGuidance("  classic:  use simple placements ");
  selDetCmd->SetGuidance("  assembly: use assembly volume "  );
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("classic");
  selDetCmd->SetCandidates("classic assembly");
  selDetCmd->AvailableForStates(PreInit,Idle);

  switchCmd = new G4UIcmdWithAString("/geom/switch",this);
  switchCmd->SetGuidance("Assign the selected geometry to G4RunManager.");
  switchCmd->SetGuidance("In case the choice is present to this command,");
  switchCmd->SetGuidance("\"/geom/select\" will be invoked and then switched.");
  switchCmd->SetParameterName("choice",true);
  switchCmd->SetDefaultValue(" ");
  switchCmd->SetCandidates("classic assembly \" \"");
  switchCmd->AvailableForStates(PreInit,Idle);

  selMatCmd = new G4UIcmdWithAString("/geom/material",this);
  selMatCmd->SetGuidance("UNUSED IN THIS VERSION...\n\n");
  selMatCmd->SetGuidance("Select Material of the SimpleBox.");
  selMatCmd->SetGuidance("  Choice : Air, Al, Pb (default)");
  selMatCmd->SetParameterName("choice",true);
  selMatCmd->SetDefaultValue("Pb");
  selMatCmd->SetCandidates("Air Al Pb");
  selMatCmd->AvailableForStates(PreInit,Idle);

  myDetector->SelectDetector(defParam="classic");
  myDetector->SelectMaterial(defParam="Pb");
}

void TstVADetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
  }
  if( command == switchCmd )
  {
    if(newValues=="classic" || newValues=="assembly")
    { myDetector->SelectDetector(newValues); }
    myDetector->SwitchDetector();
  }
  if( command == selMatCmd )
  {
    myDetector->SelectMaterial(newValues);
    //G4cout << "UNUSED IN THIS VERSION...\a" << G4endl;
  }
  return;
}

