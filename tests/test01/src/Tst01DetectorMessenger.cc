
#include "Tst01DetectorMessenger.hh"

#include "Tst01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "G4ios.hh"

Tst01DetectorMessenger::Tst01DetectorMessenger(Tst01DetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/mydet/");
  mydetDir->SetGuidance("Detector setup commands.");

  selDetCmd = new G4UIcmdWithAString("/mydet/SelectDetector",this);
  selDetCmd->SetGuidance("Select Detector Setup.");
  selDetCmd->SetGuidance("  Choice : SimpleBox / Honeycomb / Assembly");
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("SimpleBox");
  selDetCmd->SetCandidates("SimpleBox Honeycomb Assembly");
  selDetCmd->AvailableForStates(PreInit,Idle);

  switchCmd = new G4UIcmdWithAString("/mydet/SwitchDetector",this);
  switchCmd->SetGuidance("Assign the selected geometry to G4RunManager.");
  switchCmd->SetGuidance("In cese detector name is associated to this command,");
  switchCmd->
  SetGuidance("\"/mydet/SelectDetector\" will be invoked and then switched.");
  switchCmd->SetParameterName("choice",true);
  switchCmd->SetDefaultValue(" ");
  switchCmd->SetCandidates("SimpleBox Honeycomb Assembly \" \"");
  switchCmd->AvailableForStates(PreInit,Idle);

  selMatCmd = new G4UIcmdWithAString("/mydet/SelectMaterial",this);
  selMatCmd->SetGuidance("Select Material of the SimpleBox.");
  selMatCmd->SetGuidance("  Choice : Air, Al, Pb (default)");
  selMatCmd->SetParameterName("choice",true);
  selMatCmd->SetDefaultValue("Pb");
  selMatCmd->SetCandidates("Air Al Pb");
  selMatCmd->AvailableForStates(PreInit,Idle);

  // Select/Switch SCG inside Detector

  selCSGcmd = new G4UIcmdWithAString("/mydet/SelectCSG",this);
  selCSGcmd->SetGuidance("Select Detector CSG");
  selCSGcmd->SetGuidance("  Choice : Box Tubs Cons Sphere");
  selCSGcmd->SetParameterName("choice",true);
  selCSGcmd->SetDefaultValue("Box");
  selCSGcmd->SetCandidates("Box Tubs Cons Sphere");
  selCSGcmd->AvailableForStates(PreInit,Idle);

  switchCSGcmd = new G4UIcmdWithAString("/mydet/SwitchCSG",this);
  switchCSGcmd->SetGuidance("Assign the selected CSG geometry to G4RunManager.");

  switchCSGcmd->
  SetGuidance("In cese detector name is associated to this command,");

  switchCSGcmd->
  SetGuidance("\"/mydet/SelectCSG\" will be invoked and then switched.");

  switchCSGcmd->SetParameterName("choice",true);
  switchCSGcmd->SetDefaultValue(" ");
  switchCSGcmd->SetCandidates("Box Tubs Cons Sphere \" \"");
  switchCSGcmd->AvailableForStates(PreInit,Idle);

  // Select/Switch Boolean inside Detector

  selBoolCmd = new G4UIcmdWithAString("/mydet/SelectBool",this);
  selBoolCmd->SetGuidance("Select Detector Boolean");
  selBoolCmd->SetGuidance("  Choice : Intersection Union Subtraction");
  selBoolCmd->SetParameterName("choice",true);
  selBoolCmd->SetDefaultValue("Intersection");
  selBoolCmd->SetCandidates("Intersection Union Subtraction");
  selBoolCmd->AvailableForStates(PreInit,Idle);

  switchBoolCmd = new G4UIcmdWithAString("/mydet/SwitchBool",this);
  switchBoolCmd->
  SetGuidance("Assign the selected Boolean geometry to G4RunManager.");

  switchBoolCmd->
  SetGuidance("In cese detector name is associated to this command,");

  switchBoolCmd->
  SetGuidance("\"/mydet/SelectBool\" will be invoked and then switched.");

  switchBoolCmd->SetParameterName("choice",true);
  switchBoolCmd->SetDefaultValue(" ");
  switchBoolCmd->SetCandidates("Intersection Union Subtraction \" \"");
  switchBoolCmd->AvailableForStates(PreInit,Idle);

  // Default selections

  myDetector->SelectDetector(defParam="SimpleBox");
  myDetector->SelectMaterial(defParam="Pb");

  myDetector->SelectCSG(defParam="Box");
  myDetector->SelectBoolean(defParam="Intersection");
}

//////////////////////////////////////////////////////////////////////////
//
//

void Tst01DetectorMessenger::SetNewValue( G4UIcommand* command ,
                                          G4String newValues        )
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
  }
  if( command == switchCmd )
  {
    if(newValues=="SimpleBox" || newValues=="Honeycomb" || newValues=="Assembly")
    { 
      myDetector->SelectDetector(newValues); 
    }
    myDetector->SwitchDetector();
  }
  if( command == selMatCmd )
  {
    myDetector->SelectMaterial(newValues);
  }

  // Select/Switch CSG

  if( command == selCSGcmd )
  {
    myDetector->SelectCSG(newValues);
  }
  if( command == switchCSGcmd )
  {
    if( newValues=="Box"  || newValues=="Tubs" ||
        newValues=="Cons" || newValues=="Sphere"    )
    { 
      myDetector->SelectCSG(newValues); 
    }
    myDetector->SwitchCSG();
  }

  // Select/Switch Boolean

  if( command == selBoolCmd )
  {
    myDetector->SelectBoolean(newValues);
  }
  if( command == switchBoolCmd )
  {
    if( newValues=="Intersection"  || newValues=="Union"   ||
        newValues=="Subtraction"                                 )
    { 
      myDetector->SelectBoolean(newValues); 
    }
    myDetector->SwitchBoolean();
  }


  return;
}

