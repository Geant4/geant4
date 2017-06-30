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
// $Id: G4UIcontrolMessenger.cc 102561 2017-02-09 08:16:05Z gcosmo $
//

#include <stdlib.h>
#include "G4UIcontrolMessenger.hh"
#include "G4UImanager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIaliasList.hh"
#include "G4StateManager.hh"
#include "G4UIsession.hh"
#include "G4Tokenizer.hh"

#include "G4ios.hh"

G4UIcontrolMessenger::G4UIcontrolMessenger()
{
  controlDirectory = new G4UIdirectory("/control/");
  controlDirectory->SetGuidance("UI control commands.");

  macroPathCommand = new G4UIcmdWithAString("/control/macroPath",this);
  macroPathCommand->SetGuidance("Set macro search path" 
                                "with colon-separated list.");
  macroPathCommand->SetParameterName("path",false);

  ExecuteCommand = new G4UIcmdWithAString("/control/execute",this);
  ExecuteCommand->SetGuidance("Execute a macro file.");
  ExecuteCommand->SetParameterName("fileName",false);
  ExecuteCommand->SetToBeBroadcasted(false);

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
  loopCommand->SetToBeBroadcasted(false);

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
  foreachCommand->SetToBeBroadcasted(false);
  
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
  
  doublePrecCommand = new G4UIcmdWithABool("/control/useDoublePrecision",this);
  doublePrecCommand->SetGuidance("Use double precision for printing out the current parameter value(s).");
  doublePrecCommand->SetParameterName("useDoublePrecision",true);
  doublePrecCommand->SetDefaultValue(true);

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

  getValCmd = new G4UIcommand("/control/getVal",this);
  getValCmd->SetGuidance("Get the current value of the UI command and define it as an alias.");
  getValCmd->SetGuidance("Command is ignored if the UI command does not support GetCurrentValue().");
  getValCmd->SetGuidance(" Syntax : <alias_name> <UI_command> <iIdx>");
  G4UIparameter* aliName = new G4UIparameter("alias_name",'s',false);
  getValCmd->SetParameter(aliName);
  G4UIparameter* comName = new G4UIparameter("UI_command",'s',false);
  getValCmd->SetParameter(comName);
  G4UIparameter* iIdxParam = new G4UIparameter("iIdx",'i',true);
  iIdxParam->SetDefaultValue(0);
  getValCmd->SetParameter(iIdxParam);

  echoCmd = new G4UIcmdWithAString("/control/echo",this);
  echoCmd->SetGuidance("Display the aliased value.");

  shellCommand = new G4UIcmdWithAString("/control/shell",this);
  shellCommand->SetGuidance("Execute a (Unix) SHELL command.");

  ManualCommand = new G4UIcmdWithAString("/control/manual",this);
  ManualCommand->SetGuidance("Display all of sub-directories and commands.");
  ManualCommand->SetGuidance("Directory path should be given by FULL-PATH.");
  ManualCommand->SetParameterName("dirPath",true);
  ManualCommand->SetDefaultValue("/");
  ManualCommand->SetToBeBroadcasted(false);

  HTMLCommand = new G4UIcmdWithAString("/control/createHTML",this);
  HTMLCommand->SetGuidance("Generate HTML files for all of sub-directories and commands.");
  HTMLCommand->SetGuidance("Directory path should be given by FULL-PATH.");
  HTMLCommand->SetParameterName("dirPath",true);
  HTMLCommand->SetDefaultValue("/");
  HTMLCommand->SetToBeBroadcasted(false);

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
  G4UIparameter* macroFileParam = new G4UIparameter("macroFile",'s',false);
  ifCommand->SetParameter(macroFileParam);
  ifCommand->SetToBeBroadcasted(false);

  doifCommand = new G4UIcommand("/control/doif",this);
  doifCommand->SetGuidance("Execute a UI command if the expression is true.");
  doifCommand->SetGuidance(" Syntax : <double> <comp> <double> <UI_command>");
  G4UIparameter* doleftParam = new G4UIparameter("left",'d',false);
  doifCommand->SetParameter(doleftParam);
  G4UIparameter* docompParam = new G4UIparameter("comp",'s',false);
  docompParam->SetParameterCandidates("> >= < <= == !=");
  doifCommand->SetParameter(docompParam);
  G4UIparameter* dorightParam = new G4UIparameter("right",'d',false);
  doifCommand->SetParameter(dorightParam);
  G4UIparameter* comParam = new G4UIparameter("UI_command",'s',false);
  doifCommand->SetParameter(comParam);
  doifCommand->SetToBeBroadcasted(false);

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
  addCommand->SetToBeBroadcasted(false);

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
  subtractCommand->SetToBeBroadcasted(false);

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
  multiplyCommand->SetToBeBroadcasted(false);

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
  divideCommand->SetToBeBroadcasted(false);

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
  remainderCommand->SetToBeBroadcasted(false);

  strifCommand = new G4UIcommand("/control/strif",this);
  strifCommand->SetGuidance("Execute a macro file if the expression is true.");
  strifCommand->SetGuidance(" Syntax : <string> <comp> <string> <macro_file>");
  G4UIparameter* strleftParam = new G4UIparameter("left",'s',false);
  strifCommand->SetParameter(strleftParam);
  G4UIparameter* strcompParam = new G4UIparameter("comp",'s',false);
  strcompParam->SetParameterCandidates("== !=");
  strifCommand->SetParameter(strcompParam);
  G4UIparameter* strrightParam = new G4UIparameter("right",'s',false);
  strifCommand->SetParameter(strrightParam);
  G4UIparameter* strmacroFileParam = new G4UIparameter("macroFile",'s',false);
  strifCommand->SetParameter(strmacroFileParam);
  strifCommand->SetToBeBroadcasted(false);

  strdoifCommand = new G4UIcommand("/control/strdoif",this);
  strdoifCommand->SetGuidance("Execute a UI command if the expression is true.");
  strdoifCommand->SetGuidance(" Syntax : <string> <comp> <string> <UI_command>");
  G4UIparameter* strdoleftParam = new G4UIparameter("left",'s',false);
  strdoifCommand->SetParameter(strdoleftParam);
  G4UIparameter* strdocompParam = new G4UIparameter("comp",'s',false);
  strdocompParam->SetParameterCandidates("== !=");
  strdoifCommand->SetParameter(strdocompParam);
  G4UIparameter* strdorightParam = new G4UIparameter("right",'s',false);
  strdoifCommand->SetParameter(strdorightParam);
  G4UIparameter* strdomacroFileParam = new G4UIparameter("UI_command",'s',false);
  strdoifCommand->SetParameter(strdomacroFileParam);
  strdoifCommand->SetToBeBroadcasted(false);

  ifBatchCommand = new G4UIcmdWithAString("/control/ifBatch",this);
  ifBatchCommand->SetGuidance("Execute a macro file if program is running in batch mode.");
  ifBatchCommand->SetParameterName("macroFile",false);
  ifBatchCommand->SetToBeBroadcasted(false);

  ifInteractiveCommand = new G4UIcmdWithAString("/control/ifInteractive",this);
  ifInteractiveCommand->SetGuidance("Execute a macro file if program is running in interactive mode.");
  ifInteractiveCommand->SetParameterName("macroFile",false);
  ifInteractiveCommand->SetToBeBroadcasted(false);

  doifBatchCommand = new G4UIcmdWithAString("/control/doifBatch",this);
  doifBatchCommand->SetGuidance("Execute a UI command if program is running in batch mode.");
  doifBatchCommand->SetParameterName("UIcommand",false);
  doifBatchCommand->SetToBeBroadcasted(false);

  doifInteractiveCommand = new G4UIcmdWithAString("/control/doifInteractive",this);
  doifInteractiveCommand->SetGuidance("Execute a UI command if program is running in interactive mode.");
  doifInteractiveCommand->SetParameterName("UIcommand",false);
  doifInteractiveCommand->SetToBeBroadcasted(false);

}

G4UIcontrolMessenger::~G4UIcontrolMessenger()
{
  delete macroPathCommand;
  delete ExecuteCommand;
  delete suppressAbortionCommand;
  delete verboseCommand;
  delete doublePrecCommand;
  delete historyCommand;
  delete stopStoreHistoryCommand;
  delete ManualCommand;
  delete aliasCommand;
  delete unaliasCommand;
  delete listAliasCommand;
  delete getEnvCmd;
  delete getValCmd;
  delete echoCmd;
  delete shellCommand;
  delete loopCommand;
  delete foreachCommand; 
  delete HTMLCommand;
  delete maxStoredHistCommand;
  delete ifCommand;
  delete doifCommand;
  delete addCommand;
  delete subtractCommand;
  delete multiplyCommand;
  delete divideCommand;
  delete remainderCommand;
  delete strifCommand;
  delete strdoifCommand;
  delete ifBatchCommand;
  delete ifInteractiveCommand;
  delete doifBatchCommand;
  delete doifInteractiveCommand;

  delete controlDirectory;
}

void G4UIcontrolMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  G4UImanager * UI = G4UImanager::GetUIpointer();

  if( command == macroPathCommand) {
    UI-> SetMacroSearchPath(newValue);
    UI-> ParseMacroSearchPath();
  }  
  if(command==ExecuteCommand)
  {
    UI-> ExecuteMacroFile(UI-> FindMacroPath(newValue));
  }
  if(command==suppressAbortionCommand)
  {
    G4StateManager::GetStateManager()->SetSuppressAbortion(suppressAbortionCommand->GetNewIntValue(newValue));
  }
  if(command==verboseCommand)
  {
    UI->SetVerboseLevel(verboseCommand->GetNewIntValue(newValue));
  }
  if(command==doublePrecCommand)
  {
    G4UImanager::UseDoublePrecisionStr(doublePrecCommand->GetNewBoolValue(newValue));
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
  if(command==getValCmd)
  {
    G4Tokenizer next(newValue);
    G4String aliName = next();
    G4String com = next();
    G4String curVal = UI->GetCurrentValues(com);
    if(!(curVal.isNull()))
    { 
      G4String theValue = curVal;
      G4String iIdx = next();
      if(!(iIdx.isNull()))
      {
        G4int idx = StoI(iIdx);
        G4Tokenizer nextVal(curVal);
        for(G4int i = 0;i<=idx;i++)
        { theValue = nextVal(); }
      }
      G4String st = "/control/alias ";
      st += aliName + " " + theValue;
      UI->ApplyCommand(st);
    }
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
  if(command==doifCommand)
  {
    G4Tokenizer next(newValue);
    G4double l = StoD(next());
    G4String comp = next();
    G4double r = StoD(next());

    G4String c1 = next();
    G4String ca;
    while(!((ca=next()).isNull()))
    {
      c1 += " ";
      c1 += ca;
    }
    if(c1(0)=='"')
    {
      G4String strippedValue;
      if(c1(c1.length()-1)=='"')
      { strippedValue = c1(1,c1.length()-2); }
      else
      { strippedValue = c1(1,c1.length()-1); }
      c1 = strippedValue;
    }

    G4bool x = false;
    if(comp==">") x = (l>r);
    else if(comp==">=") x = (l>=r);
    else if(comp=="<") x = (l<r);
    else if(comp=="<=") x = (l<=r);
    else if(comp=="==") x = (l==r);
    else if(comp=="!=") x = (l!=r);
    if(x) UI->ApplyCommand(c1);
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
  if(command==strifCommand)
  {
    G4Tokenizer next(newValue);
    G4String l = next();
    G4String comp = next();
    G4String r = next();
    G4String mac = next();
    G4bool x = false;
    if(comp=="==") { x = (l==r); }
    else if(comp=="!=") { x = (l!=r); }
    if(x) UI->ExecuteMacroFile(mac);
  }
  if(command==strdoifCommand)
  {
    G4Tokenizer next(newValue);
    G4String l = next();
    G4String comp = next();
    G4String r = next();

    G4String c1 = next();
    G4String ca;
    while(!((ca=next()).isNull()))
    {
      c1 += " ";
      c1 += ca;
    }
    if(c1(0)=='"')
    {
      G4String strippedValue;
      if(c1(c1.length()-1)=='"')
      { strippedValue = c1(1,c1.length()-2); }
      else
      { strippedValue = c1(1,c1.length()-1); }
      c1 = strippedValue;
    }

    G4bool x = false;
    if(comp=="==") { x = (l==r); }
    else if(comp=="!=") { x = (l!=r); }
    if(x) UI->ApplyCommand(c1);
  }
  if(command==ifBatchCommand)
  {
    if(G4UIsession::InSession()==0) UI->ExecuteMacroFile(UI->FindMacroPath(newValue));
  }
  if(command==ifInteractiveCommand)
  {
    if(G4UIsession::InSession()>0) UI->ExecuteMacroFile(UI->FindMacroPath(newValue));
  }
  if(command==doifBatchCommand)
  {
    if(G4UIsession::InSession()==0) UI->ApplyCommand(newValue); 
  }
  if(command==doifInteractiveCommand)
  {
    if(G4UIsession::InSession()>0) UI->ApplyCommand(newValue); 
  }
}

G4String G4UIcontrolMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4UImanager * UI = G4UImanager::GetUIpointer();
  G4String currentValue;
  
  if( command == macroPathCommand ) {
    currentValue = UI-> GetMacroSearchPath();
  }
  if(command==verboseCommand)
  {
    currentValue = verboseCommand->ConvertToString(UI->GetVerboseLevel());
  }
  if(command==doublePrecCommand)
  {
    currentValue = doublePrecCommand->ConvertToString(G4UImanager::DoublePrecisionStr());
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


