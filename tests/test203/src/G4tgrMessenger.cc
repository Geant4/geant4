
#include "G4tgrMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

int G4tgrMessenger::theVerboseLevel = 0;

// --------------------------------------------------------------------
G4tgrMessenger::G4tgrMessenger()
{
  tgDirectory = new G4UIdirectory("/geometry/textInput/");
  tgDirectory->SetGuidance("Geometry from text file control commands.");
  verboseCmd = new G4UIcmdWithAnInteger("/geometry/textInput/verbose",this);
  verboseCmd->SetGuidance("Set Verbose level of geometry text input category.");
  verboseCmd->SetGuidance(" 0 : Silent");
  verboseCmd->SetGuidance(" 1 : info verbosity");
  verboseCmd->SetGuidance(" 2 : debug verbosity");
  verboseCmd->SetParameterName("level",false);
  verboseCmd->SetRange("level>=0");
}

// --------------------------------------------------------------------
G4tgrMessenger::~G4tgrMessenger()
{
  delete verboseCmd;
  delete tgDirectory;
}

// --------------------------------------------------------------------
void G4tgrMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == verboseCmd ){
    G4tgrMessenger::SetVerboseLevel(verboseCmd->GetNewIntValue(newValues)); 
  }
}

// --------------------------------------------------------------------
G4String G4tgrMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  if( command == verboseCmd ){ 
    cv = verboseCmd->ConvertToString(G4tgrMessenger::GetVerboseLevel()); 
  }
  return cv;
}
