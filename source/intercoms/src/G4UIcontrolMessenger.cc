//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4UIcontrolMessenger.cc,v 1.9 2002-05-15 06:51:32 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include <stdlib.h>
#include "G4UIcontrolMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIaliasList.hh"
#include "G4StateManager.hh"
#include "G4ios.hh"

G4UIcontrolMessenger::G4UIcontrolMessenger()
{
  controlDirectory = new G4UIdirectory("/control/");
  controlDirectory->SetGuidance("UI control commands.");

  ExecuteCommand = new G4UIcmdWithAString("/control/execute",this);
  ExecuteCommand->SetGuidance("Execute a macro file.");
  ExecuteCommand->SetParameterName("fileName",false);

  loopCommand = new G4UIcommand("/control/loop",this);
  loopCommand->SetGuidance("Execute a macro file more than once.");
  loopCommand->SetGuidance("Loop counter can be used as an aliased variable.");
  G4UIparameter* param1 = new G4UIparameter("macroFile",'s',false);
  loopCommand->SetParameter(param1);
  G4UIparameter* param2 = new G4UIparameter("counterName",'s',false);
  loopCommand->SetParameter(param2);
  G4UIparameter* param3 = new G4UIparameter("initialValue",'d',false);
  loopCommand->SetParameter(param3);
  G4UIparameter* param4 = new G4UIparameter("finalValue",'d',false);
  loopCommand->SetParameter(param4);
  G4UIparameter* param5 = new G4UIparameter("stepSize",'d',true);
  param5->SetDefaultValue(1.0);
  loopCommand->SetParameter(param5);

  foreachCommand = new G4UIcommand("/control/foreach",this);
  foreachCommand->SetGuidance("Execute a macro file more than once.");
  foreachCommand->SetGuidance("Loop counter can be used as an aliased variable.");
  foreachCommand->SetGuidance("Values must be separated by a space.");
  G4UIparameter* param6 = new G4UIparameter("macroFile",'s',false);
  foreachCommand->SetParameter(param6);
  G4UIparameter* param7 = new G4UIparameter("counterName",'s',false);
  foreachCommand->SetParameter(param7);
  G4UIparameter* param8 = new G4UIparameter("valueList",'s',false);
  foreachCommand->SetParameter(param8);
  
  suppressAbortionCommand = new G4UIcmdWithAnInteger("/control/suppressAbortion",this);
  suppressAbortionCommand->SetGuidance("Suppress the program abortion caused by G4Exception.");
  suppressAbortionCommand->SetGuidance("Suppression level = 0 : no suppression");
  suppressAbortionCommand->SetGuidance("                  = 1 : suppress during EventProc state");
  suppressAbortionCommand->SetGuidance("                  = 2 : full suppression, i.e. no abortion by G4Exception");
  suppressAbortionCommand->SetGuidance("When abortion is suppressed, you will get error messages issued by G4Exception,");
  suppressAbortionCommand->SetGuidance("and there is NO guarantee for the correct result after the G4Exception error message.");
  suppressAbortionCommand->SetParameterName("level",true);
  suppressAbortionCommand->SetRange("level >= 0 && level <= 2");
  suppressAbortionCommand->SetDefaultValue(0);

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

  aliasCommand = new G4UIcommand("/control/alias",this);
  aliasCommand->SetGuidance("Set an alias.");
  aliasCommand->SetGuidance("String can be aliased by this command.");
  aliasCommand->SetGuidance("The string may contain one or more spaces,");
  aliasCommand->SetGuidance("the string must be enclosed by double quotes (\").");
  aliasCommand->SetGuidance("To use an alias, enclose the alias name with");
  aliasCommand->SetGuidance("parenthis \"{\" and \"}\".");
  G4UIparameter* aliasNameParam = new G4UIparameter("aliasName",'s',false);
  aliasCommand->SetParameter(aliasNameParam);
  G4UIparameter* aliasValueParam = new G4UIparameter("aliasValue",'s',false);
  aliasCommand->SetParameter(aliasValueParam);

  unaliasCommand = new G4UIcmdWithAString("/control/unalias",this);
  unaliasCommand->SetGuidance("Remove an alias.");
  unaliasCommand->SetParameterName("aliasName",false);

  listAliasCommand = new G4UIcmdWithoutParameter("/control/listAlias",this);
  listAliasCommand->SetGuidance("List aliases.");

  shellCommand = new G4UIcmdWithAString("/control/shell",this);
  shellCommand->SetGuidance("Execute a (Unix) SHELL command.");

  ManualCommand = new G4UIcmdWithAString("/control/manual",this);
  ManualCommand->SetGuidance("Display all of sub-directories and commands.");
  ManualCommand->SetGuidance("Directory path should be given by FULL-PATH.");
  ManualCommand->SetParameterName("dirPath",true);
  ManualCommand->SetDefaultValue("/");

  HTMLCommand = new G4UIcmdWithAString("/control/createHTML",this);
  HTMLCommand->SetGuidance("Generate HTML files for all of sub-directories and commands.");
  HTMLCommand->SetGuidance("Directory path should be given by FULL-PATH.");
  HTMLCommand->SetParameterName("dirPath",true);
  HTMLCommand->SetDefaultValue("/");

  maxStoredHistCommand = new G4UIcmdWithAnInteger("/control/maximumStoredHistory",this);
  maxStoredHistCommand->SetGuidance("Set maximum number of stored UI commands.");
  maxStoredHistCommand->SetParameterName("max",true);
  maxStoredHistCommand->SetDefaultValue(20);
}

G4UIcontrolMessenger::~G4UIcontrolMessenger()
{
  delete ExecuteCommand;
  delete suppressAbortionCommand;
  delete verboseCommand;
  delete historyCommand;
  delete stopStoreHistoryCommand;
  delete ManualCommand;
  delete aliasCommand;
  delete unaliasCommand;
  delete listAliasCommand;
  delete shellCommand;
  delete loopCommand;
  delete foreachCommand; 
  delete HTMLCommand;
  delete maxStoredHistCommand;
  delete controlDirectory;
}

void G4UIcontrolMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4UImanager * UI = G4UImanager::GetUIpointer();
  
  if(command==ExecuteCommand)
  {
    UI->ExecuteMacroFile(newValue);
  }
  if(command==suppressAbortionCommand)
  {
    G4StateManager::GetStateManager()->SetSuppressAbortion(suppressAbortionCommand->GetNewIntValue(newValue));
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
  if(command==aliasCommand)
  {
    UI->SetAlias(newValue);
  }
  if(command==unaliasCommand)
  {
    UI->RemoveAlias(newValue);
  }
  if(command==listAliasCommand)
  {
    UI->ListAlias();
  }
  if(command==shellCommand)
  {
    system(newValue);
  }
  if(command==loopCommand)
  {
    UI->LoopS(newValue);
  }
  if(command==foreachCommand)
  {
    UI->ForeachS(newValue);
  }
  if(command==HTMLCommand)
  {
    UI->CreateHTML(newValue);
  }
  if(command==maxStoredHistCommand)
  {
    UI->SetMaxHistSize(maxStoredHistCommand->GetNewIntValue(newValue));
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
  if(command==suppressAbortionCommand)
  {
    currentValue = suppressAbortionCommand->ConvertToString(G4StateManager::GetStateManager()->GetSuppressAbortion());
  }
  if(command==maxStoredHistCommand)
  {
    currentValue = maxStoredHistCommand->ConvertToString(UI->GetMaxHistSize());
  }
  
  return currentValue;
}


