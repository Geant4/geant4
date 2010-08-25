//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4UIcontrolMessenger.cc,v 1.12 2010-08-25 06:09:57 asaim Exp $
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
#include "G4Tokenizer.hh"

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

  getEnvCmd = new G4UIcmdWithAString("/control/getEnv",this);
  getEnvCmd->SetGuidance("Get a shell environment variable and define it as an alias.");

  echoCmd = new G4UIcmdWithAString("/control/echo",this);
  echoCmd->SetGuidance("Display the aliased value.");

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

  ifCommand = new G4UIcommand("/control/if",this);
  ifCommand->SetGuidance("Execute a macro file if the expression is true.");
  ifCommand->SetGuidance(" Syntax : <double> <comp> <double> <macro_file>");
  G4UIparameter* leftParam = new G4UIparameter("left",'d',false);
  ifCommand->SetParameter(leftParam);
  G4UIparameter* compParam = new G4UIparameter("comp",'s',false);
  compParam->SetParameterCandidates("> >= < <= == !=");
  ifCommand->SetParameter(compParam);
  G4UIparameter* rightParam = new G4UIparameter("right",'d',false);
  ifCommand->SetParameter(rightParam);
  G4UIparameter* macroFileParam = new G4UIparameter("aliasValue",'s',false);
  ifCommand->SetParameter(macroFileParam);

  addCommand = new G4UIcommand("/control/add",this);
  addCommand->SetGuidance("Define a new alias as the sum of two values.");
  addCommand->SetGuidance(" Syntax : <new_alias> <value1> <value2>");
  addCommand->SetGuidance(" <new_alias> may be an already existing alias. If it is the case,");
  addCommand->SetGuidance(" aliased value is alternated.");
  G4UIparameter* newAlias1 = new G4UIparameter("new_alias",'s',false);
  addCommand->SetParameter(newAlias1);
  G4UIparameter* val1a = new G4UIparameter("value1",'d',false);
  addCommand->SetParameter(val1a);
  G4UIparameter* val1b = new G4UIparameter("value2",'d',false);
  addCommand->SetParameter(val1b);

  subtractCommand = new G4UIcommand("/control/subtract",this);
  subtractCommand->SetGuidance("Define a new alias as the subtraction of two values.");
  subtractCommand->SetGuidance(" Syntax : <new_alias> <value1> <value2>");
  subtractCommand->SetGuidance(" <new_alias> may be an already existing alias. If it is the case,");
  subtractCommand->SetGuidance(" aliased value is alternated.");
  G4UIparameter* newAlias2 = new G4UIparameter("new_alias",'s',false);
  subtractCommand->SetParameter(newAlias2);
  G4UIparameter* val2a = new G4UIparameter("value1",'d',false);
  subtractCommand->SetParameter(val2a);
  G4UIparameter* val2b = new G4UIparameter("value2",'d',false);
  subtractCommand->SetParameter(val2b);

  multiplyCommand = new G4UIcommand("/control/multiply",this);
  multiplyCommand->SetGuidance("Define a new alias as the multiplification of two values.");
  multiplyCommand->SetGuidance(" Syntax : <new_alias> <value1> <value2>");
  multiplyCommand->SetGuidance(" <new_alias> may be an already existing alias. If it is the case,");
  multiplyCommand->SetGuidance(" aliased value is alternated.");
  G4UIparameter* newAlias3 = new G4UIparameter("new_alias",'s',false);
  multiplyCommand->SetParameter(newAlias3);
  G4UIparameter* val3a = new G4UIparameter("value1",'d',false);
  multiplyCommand->SetParameter(val3a);
  G4UIparameter* val3b = new G4UIparameter("value2",'d',false);
  multiplyCommand->SetParameter(val3b);

  divideCommand = new G4UIcommand("/control/divide",this);
  divideCommand->SetGuidance("Define a new alias as the division of two values.");
  divideCommand->SetGuidance(" Syntax : <new_alias> <value1> <value2>");
  divideCommand->SetGuidance(" <new_alias> may be an already existing alias. If it is the case,");
  divideCommand->SetGuidance(" aliased value is alternated.");
  G4UIparameter* newAlias4 = new G4UIparameter("new_alias",'s',false);
  divideCommand->SetParameter(newAlias4);
  G4UIparameter* val4a = new G4UIparameter("value1",'d',false);
  divideCommand->SetParameter(val4a);
  G4UIparameter* val4b = new G4UIparameter("value2",'d',false);
  val4b->SetParameterRange("value2 != 0.");
  divideCommand->SetParameter(val4b);

  remainderCommand = new G4UIcommand("/control/remainder",this);
  remainderCommand->SetGuidance("Define a new alias as the remainder of two values.");
  remainderCommand->SetGuidance(" Syntax : <new_alias> <value1> <value2>");
  remainderCommand->SetGuidance(" <new_alias> may be an already existing alias. If it is the case,");
  remainderCommand->SetGuidance(" aliased value is alternated.");
  G4UIparameter* newAlias5 = new G4UIparameter("new_alias",'s',false);
  remainderCommand->SetParameter(newAlias5);
  G4UIparameter* val5a = new G4UIparameter("value1",'i',false);
  remainderCommand->SetParameter(val5a);
  G4UIparameter* val5b = new G4UIparameter("value2",'i',false);
  val4b->SetParameterRange("value2 != 0");
  remainderCommand->SetParameter(val5b);

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
  delete getEnvCmd;
  delete echoCmd;
  delete shellCommand;
  delete loopCommand;
  delete foreachCommand; 
  delete HTMLCommand;
  delete maxStoredHistCommand;
  delete ifCommand;
  delete addCommand;
  delete subtractCommand;
  delete multiplyCommand;
  delete divideCommand;
  delete remainderCommand;

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
  if(command==getEnvCmd)
  {
    if(getenv(newValue))
    { 
      G4String st = "/control/alias ";
      st += newValue;
      st += " ";
      st += getenv(newValue);
      UI->ApplyCommand(st);
    }
    else
    { G4cerr << "<" << newValue << "> is not defined as a shell variable. Command ignored." << G4endl; }
  }
  if(command==echoCmd)
  { G4cout << UI->SolveAlias(newValue) << G4endl; }
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
  if(command==ifCommand)
  {
    G4Tokenizer next(newValue);
    G4double l = StoD(next());
    G4String comp = next();
    G4double r = StoD(next());
    G4String mac = next();
    G4bool x = false;
    if(comp==">") x = (l>r);
    else if(comp==">=") x = (l>=r);
    else if(comp=="<") x = (l<r);
    else if(comp=="<=") x = (l<=r);
    else if(comp=="==") x = (l==r);
    else if(comp=="!=") x = (l!=r);
    if(x) UI->ExecuteMacroFile(mac);
  }
  if(command==addCommand)
  {
    G4Tokenizer next(newValue);
    G4String newA = next();
    G4double l = StoD(next());
    G4double r = StoD(next());
    G4String st = "/control/alias ";
    st += newA;
    st += " ";
    st += DtoS(l+r);
    UI->ApplyCommand(st);
  }
  if(command==subtractCommand)
  {
    G4Tokenizer next(newValue);
    G4String newA = next();
    G4double l = StoD(next());
    G4double r = StoD(next());
    G4String st = "/control/alias ";
    st += newA;
    st += " ";
    st += DtoS(l-r);
    UI->ApplyCommand(st);
  }
  if(command==multiplyCommand)
  {
    G4Tokenizer next(newValue);
    G4String newA = next();
    G4double l = StoD(next());
    G4double r = StoD(next());
    G4String st = "/control/alias ";
    st += newA;
    st += " ";
    st += DtoS(l*r);
    UI->ApplyCommand(st);
  }
  if(command==divideCommand)
  {
    G4Tokenizer next(newValue);
    G4String newA = next();
    G4double l = StoD(next());
    G4double r = StoD(next());
    G4String st = "/control/alias ";
    st += newA;
    st += " ";
    st += DtoS(l/r);
    UI->ApplyCommand(st);
  }
  if(command==remainderCommand)
  {
    G4Tokenizer next(newValue);
    G4String newA = next();
    G4int l = StoI(next());
    G4int r = StoI(next());
    G4String st = "/control/alias ";
    st += newA;
    st += " ";
    st += DtoS(l%r);
    UI->ApplyCommand(st);
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


