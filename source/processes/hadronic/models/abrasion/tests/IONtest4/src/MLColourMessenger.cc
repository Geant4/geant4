////////////////////////////////////////////////////////////////////////////////
//
#include "MLColourMessenger.hh"

#include "MLColour.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithoutParameter.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLColourMessenger::MLColourMessenger(MLColour * MLCol)
:colourManager(MLCol)
{ 
  ColourDir = new G4UIdirectory("/geometry/colour/");
  ColourDir->SetGuidance("Controls for defining geometry colours." );

  AddCmd = new G4UIcommand("/geometry/colour/add",this);
  AddCmd->SetGuidance("Add a colour by its name and RGB intensities.");
  AddCmd->SetGuidance("  0. <= RGBs <= 1. ");
  //
  G4UIparameter* ColName = new G4UIparameter("colour",'s',false);
  ColName->SetGuidance("colour name");
  AddCmd->SetParameter(ColName);
  //
  G4UIparameter* Red = new G4UIparameter("red",'d',false);
  Red->SetGuidance("red component");
  Red->SetParameterRange("red >= 0. &&  red <= 1.");
  AddCmd->SetParameter(Red);
 //
  G4UIparameter* Green = new G4UIparameter("green",'d',false);
  Green->SetGuidance("green component");
  Green->SetParameterRange("green >= 0. && green <= 1.");
  AddCmd->SetParameter(Green);
 //
  G4UIparameter* Blue = new G4UIparameter("blue",'d',false);
  Blue->SetGuidance("blue component");
  Blue->SetParameterRange("blue >= 0. && blue <= 1.");
  AddCmd->SetParameter(Blue);

  AddCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ListCmd = new G4UIcmdWithoutParameter("/geometry/colour/list",this);
  ListCmd->SetGuidance("List the colours defined.");
  ListCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}
////////////////////////////////////////////////////////////////////////////////
//
MLColourMessenger::~MLColourMessenger()
{
  delete ColourDir;
  delete AddCmd;
  delete ListCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void MLColourMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{    
              
  if (command == ListCmd) {
    colourManager->ListColour();
   
  } else if (command == AddCmd) {
    G4double r,g,b;
    char mat[30];
    const char* t = newValue;
    std::istrstream is((char*)t);
    is >>mat >>r >>g >>b;
    G4String name=mat;
    colourManager->AddColour(name,r,g,b);

  }
}
////////////////////////////////////////////////////////////////////////////////
