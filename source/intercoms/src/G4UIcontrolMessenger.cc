// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIcontrolMessenger.cc,v 1.1 1999-01-07 16:09:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIcontrolMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4ios.hh"

G4UIcontrolMessenger::G4UIcontrolMessenger()
{
  controlDirectory = new G4UIdirectory("/control/");
  controlDirectory->SetGuidance("UI control commands.");

  ExecuteCommand = new G4UIcmdWithAString("/control/execute",this);
  ExecuteCommand->SetGuidance("Execute a macro file.");
  ExecuteCommand->SetParameterName("fileName",false);
  
  verboseCommand = new G4UIcmdWithAnInteger("/control/verbose",this);
  verboseCommand->SetGuidance("Applied command will also be shown on screen.");
  verboseCommand->SetGuidance("This command is useful with MACRO file.");
  verboseCommand->SetGuidance("  0 : silent");
  verboseCommand->SetGuidance("  1 : only the valid commands are shown.");
  verboseCommand->SetGuidance("  2 : comment lines are also shown (default).");
  verboseCommand->SetParameterName("switch",true);
  verboseCommand->SetRange("switch >= 0 && switch <=2");
  verboseCommand->SetDefaultValue(2);
  
  historyCommand = new G4UIcmdWithAString("/control/saveHistory",this);
  historyCommand->SetGuidance("Store command history to a file.");
  historyCommand->SetGuidance("Defaul file name is G4history.macro.");
  historyCommand->SetParameterName("fileName",true);
  historyCommand->SetDefaultValue("G4History.macro");
  
  stopStoreHistoryCommand 
    = new G4UIcmdWithoutParameter("/control/stopSavingHistory",this);
  stopStoreHistoryCommand->SetGuidance("Stop saving history file.");

  ManualCommand = new G4UIcmdWithAString("/control/manual",this);
  ManualCommand->SetGuidance("Display all of sub-directories and commands.");
  ManualCommand->SetGuidance("Directory path should be given by FULL-PATH.");
  ManualCommand->SetParameterName("dirPath",true);
  ManualCommand->SetDefaultValue("/");
}

G4UIcontrolMessenger::~G4UIcontrolMessenger()
{
  delete ExecuteCommand;
  delete verboseCommand;
  delete historyCommand;
  delete stopStoreHistoryCommand;
  delete ManualCommand;
  
  delete controlDirectory;
}

void G4UIcontrolMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4UImanager * UI = G4UImanager::GetUIpointer();
  
  if(command==ExecuteCommand)
  {
    UI->ExecuteMacroFile(newValue);
  }
  if(command==verboseCommand)
  {
    UI->SetVerboseLevel(verboseCommand->GetNewIntValue(newValue));
  }
  if(command==historyCommand)
  {
  	UI->StoreHistory(newValue);
  }
  if(command==stopStoreHistoryCommand)
  {
  	UI->StoreHistory(false);
  }
  if(command==ManualCommand)
  {
    UI->ListCommands(newValue);
  }
}

G4String G4UIcontrolMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4UImanager * UI = G4UImanager::GetUIpointer();
  G4String currentValue;
  
  if(command==verboseCommand)
  {
    currentValue = verboseCommand->ConvertToString(UI->GetVerboseLevel());
  }
  
  return currentValue;
}


