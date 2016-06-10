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
// $Id: G4UIterminal.cc 81292 2014-05-26 09:32:37Z gcosmo $
//
// ====================================================================
//   G4UIterminal.cc
//
// ====================================================================
#include "G4Types.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4UIterminal.hh"
#include "G4UIcsh.hh"
#include <sstream>

#ifndef WIN32
#include <signal.h>
#endif

// ====================================================================
// signal handler for soft-abort
// ====================================================================

static G4ThreadLocal G4VUIshell* theshell= 0;

#ifndef WIN32

extern "C" {

////////////////////////////////
static void SignalHandler(G4int)
////////////////////////////////
{
  G4StateManager* stateManager= G4StateManager::GetStateManager();
  G4ApplicationState state= stateManager-> GetCurrentState();

  if(state==G4State_GeomClosed || state==G4State_EventProc) {
    G4cout << "aborting Run ...";
    G4UImanager::GetUIpointer()->ApplyCommand("/run/abort");
    G4cout << G4endl;
  } else {
    G4cout << G4endl
           << "Session terminated." << G4endl;
    theshell-> ResetTerminal();
    G4Exception("G4UIterminal::SignalHandler()",
                "UI0001",
                FatalException, 
                "KeyboardInterrput with Ctrl-C");
  }

  // for original Unix / System V
  signal(SIGINT, SignalHandler);
}

}
#endif

// ====================================================================
//
// class description
//
// ====================================================================

///////////////////////////////////////////////////////////
G4UIterminal::G4UIterminal(G4VUIshell* aShell, G4bool qsig)
///////////////////////////////////////////////////////////
{
  UI= G4UImanager::GetUIpointer();
  UI-> SetSession(this);
  UI-> SetCoutDestination(this);

  iExit= FALSE;
  iCont= FALSE;

  if(aShell) shell= aShell;
  else shell= new G4UIcsh;
  theshell= shell; // locally stored for the signal handler

  // add signal handler
  if(qsig) {
#ifndef WIN32
  signal(SIGINT, SignalHandler);
#endif
  }
}

/////////////////////////////
G4UIterminal::~G4UIterminal() 
/////////////////////////////
{ 
  if(shell) delete shell;

  if(G4UImanager::GetUIpointer()) {
    UI-> SetSession(NULL);
    UI-> SetCoutDestination(NULL);
  }
}


////////////////////////////////////////////////////
void G4UIterminal::SetPrompt(const G4String& prompt) 
////////////////////////////////////////////////////
{
  shell-> SetPrompt(prompt);
}

/////////////////////////////////////////
G4UIsession* G4UIterminal::SessionStart()
/////////////////////////////////////////
{
  iExit= TRUE;

  G4String newCommand= GetCommand();
  while(iExit){
    ExecuteCommand(newCommand);
    newCommand= GetCommand();
  }
  return NULL;
}

/////////////////////////////////////////////////////////
void G4UIterminal::PauseSessionStart(const G4String& msg)
/////////////////////////////////////////////////////////
{
  iCont= TRUE;

  G4String newCommand= GetCommand(msg);
  while(iCont){
    ExecuteCommand(newCommand);
    newCommand= GetCommand(msg);
  }
}

///////////////////////////////////////////////////////////
void G4UIterminal::ExecuteCommand(const G4String& aCommand)
///////////////////////////////////////////////////////////
{
  if(aCommand.length()<2) return;

  G4int returnVal = UI-> ApplyCommand(aCommand);

  G4int paramIndex = returnVal % 100;
  // 0 - 98 : paramIndex-th parameter is invalid
  // 99     : convination of parameters is invalid
  G4int commandStatus = returnVal - paramIndex;

  G4UIcommand* cmd = 0;
  if(commandStatus!=fCommandSucceeded)
  { cmd = FindCommand(aCommand); }

  switch(commandStatus) {
  case fCommandSucceeded:
    break;
  case fCommandNotFound:
    G4cerr << "command <" << UI->SolveAlias(aCommand) << "> not found" << G4endl;
    if( aCommand.index("@@") != G4String::npos) {
      G4cout << "@@G4UIterminal" << G4endl;
    }
    break;
  case fIllegalApplicationState:
    G4cerr << "illegal application state -- command refused" << G4endl;
    break;
  case fParameterOutOfRange:
    // if(paramIndex<99) {
    //   G4cerr << "Parameter is out of range (index " << paramIndex << ")" << G4endl;
    //   G4cerr << "Allowed range : " << cmd->GetParameter(paramIndex)->GetParameterRange() << G4endl;
    // } else {
    //   G4cerr << "Parameter is out of range" << G4endl;
    //   G4cerr << "Allowed range : " << cmd->GetRange() << G4endl;
    // }
    break;
  case fParameterOutOfCandidates:
    G4cerr << "Parameter is out of candidate list (index " << paramIndex << ")" << G4endl;
    G4cerr << "Candidates : " << cmd->GetParameter(paramIndex)->GetParameterCandidates() << G4endl;
    break;
  case fParameterUnreadable:
    G4cerr << "Parameter is wrong type and/or is not omittable (index " << paramIndex << ")" << G4endl;
    break;
  case fAliasNotFound:
  default:
    G4cerr << "command refused (" << commandStatus << ")" << G4endl;
  }
}

//////////////////////////////////////////////////
G4String G4UIterminal::GetCommand(const char* msg)
//////////////////////////////////////////////////
{
  G4String newCommand;
  G4String nullString;

  newCommand= shell-> GetCommandLineString(msg);

  G4String nC = newCommand.strip(G4String::leading);
  if( nC.length() == 0 ) {
    newCommand= nullString;

  } else if( nC(0) == '#' ) {  
    G4cout << nC << G4endl;
    newCommand= nullString;

  } else if(nC=="ls" || nC.substr(0,3)=="ls " ) {  // list commands
    ListDirectory(nC); 
    newCommand= nullString;

  } else if(nC=="lc" || nC.substr(0,3)=="lc " ) {  // ... by shell
    shell-> ListCommand(nC.remove(0,2)); 
    newCommand= nullString;

  } else if(nC == "pwd") { // show current directory
    G4cout << "Current Command Directory : " 
	   << GetCurrentWorkingDirectory() << G4endl; 
    newCommand= nullString;

  } else if(nC == "cwd") { // ... by shell
    shell-> ShowCurrentDirectory();
    newCommand= nullString;

  } else if(nC == "cd" || nC.substr(0,3) == "cd ") {  // "cd"
    ChangeDirectoryCommand(nC); 
    shell-> SetCurrentDirectory(GetCurrentWorkingDirectory());
    newCommand= nullString;

  } else if(nC == "help" || nC.substr(0,5) == "help ") {  // "help"
    TerminalHelp(nC);
    newCommand= nullString;

  } else if(nC(0) == '?') {   // "show current value of a parameter"
    ShowCurrent(nC);
    newCommand= nullString;

  } else if(nC == "hist" || nC == "history") {     // "hist/history"
    G4int nh= UI-> GetNumberOfHistory();
    for (G4int i=0; i<nh; i++) { 
      G4cout << i << ": " << UI->GetPreviousCommand(i) << G4endl; 
    }
    newCommand= nullString;

  } else if(nC(0) == '!') {   // "!"
    G4String ss= nC(1, nC.length()-1);
    G4int vl;
    const char* tt= ss;
    std::istringstream is(tt);
    is >> vl;
    G4int nh= UI-> GetNumberOfHistory();
    if(vl>=0 && vl<nh) { 
      newCommand= UI-> GetPreviousCommand(vl); 
      G4cout << newCommand << G4endl;
    } else { 
      G4cerr << "history " << vl << " is not found." << G4endl; 
      newCommand= nullString;
    }

  } else if(nC == "exit") {   // "exit"
    if(iCont) { 
      G4cout << "You are now processing RUN." << G4endl;
      G4cout << "Please abort it using \"/run/abort\" command first" << G4endl;
      G4cout << " and use \"continue\" command until the application" 
	     << G4endl;
      G4cout << " becomes to Idle." << G4endl;
    } else {
      iExit= FALSE;
      newCommand= nullString;
    }

  } else if( nC == "cont" || nC == "continue"){     // "cont/continu"
    iCont= FALSE;
    newCommand= nullString;

  } else if( nC.empty() ){ // NULL command
    newCommand= nullString;
    
  } else {
  }

  return ModifyToFullPathCommand(newCommand);
}


/////////////////////////////////////////////////////////////
G4int G4UIterminal::ReceiveG4cout(const G4String& coutString)
/////////////////////////////////////////////////////////////
{
  std::cout << coutString << std::flush;
  return 0;
}

/////////////////////////////////////////////////////////////
G4int G4UIterminal::ReceiveG4cerr(const G4String& cerrString)
/////////////////////////////////////////////////////////////
{
  std::cerr << cerrString << std::flush;
  return 0;
}

///////////////////////////////////////////////
G4bool G4UIterminal::GetHelpChoice(G4int& aInt)
///////////////////////////////////////////////
{
  G4cin >> aInt;
  if(!G4cin.good()){
    G4cin.clear();
    G4cin.ignore(30,'\n');
    return FALSE;
  }
  return TRUE;
}

///////////////////////////////////
void G4UIterminal::ExitHelp() const
///////////////////////////////////
{
  char temp[100];
  G4cin.getline(temp, 100);
}

