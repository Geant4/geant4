// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIterminal.cc,v 1.1 1999-01-07 16:09:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4UIterminal.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4ios.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4UIterminal::G4UIterminal()
{
  UI = G4UImanager::GetUIpointer();
  UI->SetSession(this);
  UI->SetCoutDestination(this);
  G4StateManager * statM = G4StateManager::GetStateManager();
  promptCharacter = statM->GetStateString(statM->GetCurrentState())
                    + G4String("> ");
  iExit = false;
  iCont = false;
}

G4UIterminal::~G4UIterminal() 
{ 
 if( G4UImanager::GetUIpointer() != 0)
  {
    UI->SetSession(NULL);
    UI->SetCoutDestination(NULL);
    //    G4cout << "Terminal session deleted" << endl;
  }
}

G4UIsession * G4UIterminal::SessionStart()
{
  iExit = true;
  G4StateManager * statM = G4StateManager::GetStateManager();
  promptCharacter = statM->GetStateString(statM->GetCurrentState())
                  + G4String("> ");
  G4String newCommand = GetCommand();
  while( iExit )
  {
    ExecuteCommand(newCommand);
    promptCharacter = statM->GetStateString(statM->GetCurrentState())
                    + G4String("> ");
    newCommand = GetCommand();
  }
  return NULL;
}

void G4UIterminal::PauseSessionStart(G4String msg)
{
  promptCharacter = msg;
  iCont = true;
  G4String newCommand = GetCommand();
  while( iCont )
  {
    ExecuteCommand(newCommand);
    newCommand = GetCommand();
  }
}

G4int G4UIterminal::ReceiveG4cout(G4String coutString)
{
  cout << coutString << flush;
  return 0;
}

G4int G4UIterminal::ReceiveG4cerr(G4String cerrString)
{
  cerr << cerrString << flush;
  return 0;
}

void G4UIterminal::ExecuteCommand(G4String aCommand)
{
  if(aCommand.length()<2) return;
  int commandStatus = UI->ApplyCommand(aCommand);
  switch(commandStatus)
  {
    case fCommandSucceeded:
      break;
    case fCommandNotFound:
      G4cerr << "command <" << aCommand << "> not found" << endl;
       if ( aCommand.index("@@") != RW_NPOS) {
       G4cout << "@@G4UIterminal" << endl;
       }
      break;
    case fIllegalApplicationState:
      G4cerr << "illegal application state -- command refused" << endl;
      break;
    case fParameterOutOfRange:
    case fParameterUnreadable:
    case fParameterOutOfCandidates:
    default:
      G4cerr << "command refused (" << commandStatus << ")" << endl;
  }
}

G4String G4UIterminal::GetCommand()
{
  G4String newCommand;
  G4String nullString;
  while( 1 )
  {
    G4cout << promptCharacter;
    G4cout.flush();
    newCommand.readLine( cin );
    if (!cin.good()) { cin.clear(); newCommand = nullString; iExit=false; break; }
    while( newCommand(newCommand.length()-1) == '_' )
    {
      G4String newLine;
      newCommand.remove(newCommand.length()-1);
      newLine.readLine( cin );
      if (!cin.good()) { cin.clear(); newCommand = nullString; iExit=false; break; }
      newCommand.append(newLine);
    }

    G4String nC = newCommand.strip(G4String::leading);
    // -------------------- nC.toUpper();
    if( nC(0) == '#' )
    { G4cout << nC << endl; }
    else if( nC == "ls" || nC(0,3) == "ls " )
    { ListDirectory( nC ); }
    else if( nC == "pwd" )
    { 
      G4cout << "Current Working Directory : " 
           << GetCurrentWorkingDirectory() << endl; 
    }
    else if( nC == "cd" || nC(0,3) == "cd " )
    { ChangeDirectoryCommand( nC ); }
    else if( nC == "help" || nC(0,5) == "help " )
    { TerminalHelp( nC ); }
    else if( nC(0) == '?' )
    { ShowCurrent( nC ); }
    else if( nC == "hist" || nC == "history" )
    {
      G4int nh = UI->GetNumberOfHistory();
      for(int i=0;i<nh;i++)
      { G4cout << i << ": " << UI->GetPreviousCommand(i) << endl; }
    }
    else if( nC(0) == '!' )
    {
      G4String ss = nC(1,nC.length()-1);
      G4int vl;
      const char* tt = ss;
      istrstream is((char*)tt);
      is >> vl;
      G4int nh = UI->GetNumberOfHistory();
      if(vl>=0 && vl<nh)
      { 
        newCommand = UI->GetPreviousCommand(vl); 
        G4cout << newCommand << endl;
        break;
      }
      else
      { G4cerr << "history " << vl << " is not found." << endl; }
    }
    else if( nC == "exit" )
    { 
      if( iCont )
      { 
        G4cout << "You are now processing RUN." << endl;
        G4cout << "Please abort it using \"/run/abort\" command first" << endl;
        G4cout << " and use \"continue\" command until the application" << endl;
        G4cout << " becomes to Idle." << endl;
      }
      else
      {
        iExit = false;
        newCommand = nullString;
        break;
      }
    }
    else if( nC == "cont" || nC == "continue" )
    { 
      iCont = false;
      newCommand = nullString;
      break;
    }
    else
    { break; }
  }
  return ModifyToFullPathCommand(newCommand);
}

void G4UIterminal::ShowCurrent( G4String newCommand )
{
  G4String comString = newCommand(1,newCommand.length()-1);
  G4String theCommand = ModifyToFullPathCommand(comString);
  G4String curV = UI->GetCurrentValues(theCommand);
  if( ! curV.isNull() )
  { G4cout << "Current value(s) of the parameter(s) : " << curV << endl; }
}

void G4UIterminal::ChangeDirectoryCommand( G4String newCommand )
{
  G4String prefix;
  if( newCommand.length() <= 3 )
  { prefix = "/"; }
  else
  {
    G4String aNewPrefix = newCommand(3,newCommand.length()-3);
    prefix = aNewPrefix.strip(G4String::both);
  }
  if(!ChangeDirectory(prefix))
  { G4cout << "directory <" << prefix << "> not found." << endl; }
}

void G4UIterminal::ListDirectory( G4String newCommand )
{
  G4String targetDir;
  if( newCommand.length() <= 3 )
  { targetDir = GetCurrentWorkingDirectory(); }
  else
  {
    G4String newPrefix = newCommand(3,newCommand.length()-3);
    targetDir = newPrefix.strip(G4String::both);
  }
  G4UIcommandTree* commandTree 
    = FindDirectory(ModifyToFullPathCommand(targetDir));
  if( commandTree == NULL )
  { G4cout << "Directory <" << targetDir << "> is not found." << endl; }
  else
  { commandTree->ListCurrent(); }
}

void G4UIterminal::TerminalHelp(G4String newCommand)
{
  G4UIcommandTree * treeTop = UI->GetTree();
  int i = newCommand.index(" ");
  if( i != RW_NPOS )
  {
    G4String newValue = newCommand(i+1,newCommand.length()-(i+1));
    newValue.strip(G4String::both);
    G4String targetCom = ModifyToFullPathCommand( newValue );
    G4UIcommand* theCommand = treeTop->FindPath( targetCom );
    if( theCommand != NULL ) 
    { 
      theCommand->List();
      return;
    }
    else
    {
      G4cout << "Command <" << newValue << " is not found." << endl;
      return;
    }
  }

  G4UIcommandTree * floor[10];
  floor[0] = treeTop;
  int iFloor = 0;
  int prefixIndex = 1;
  G4String prefix = GetCurrentWorkingDirectory();
  while( prefixIndex < prefix.length()-1 )
  {
    int ii = prefix.index("/",prefixIndex);
    floor[iFloor+1] = 
      floor[iFloor]->GetTree(prefix(0,ii+1));
    prefixIndex = ii+1;
    iFloor++;
  }
  floor[iFloor]->ListCurrentWithNum();
  // 1998 Oct 2 non-number input
  while(1){
    G4cout << endl << "Type the number ( 0:end, -n:n level back ) : "<<flush;
    cin >> i;
    if(!cin.good()){
      cin.clear();
      cin.ignore(30,'\n');
      G4cout << endl << "Not a number, once more" << endl; continue;}
    else if( i < 0 ){
      iFloor += i;
      if( iFloor < 0 ) iFloor = 0;
      floor[iFloor]->ListCurrentWithNum(); continue;}
    else if(i == 0) { break;}
    else if( i > 0 ) {
      int n_tree = floor[iFloor]->GetTreeEntry();
      if( i > n_tree )
      { 
        if( i <= n_tree + floor[iFloor]->GetCommandEntry() )
        { 
          floor[iFloor]->GetCommand(i-n_tree)->List(); 
          //iFloor++;
        }
      }
      else
      {
        floor[iFloor+1] = floor[iFloor]->GetTree(i);
        iFloor++;
        floor[iFloor]->ListCurrentWithNum();
      }
    }
    //G4cout << endl << "Type the number ( 0:end, -n:n level back ) : ";

  }
  G4cout << "Exit from HELP." << endl << endl;
  G4cout << endl;
  // cin.flush();
  char temp[100];
  cin.getline( temp, 100 );
}


