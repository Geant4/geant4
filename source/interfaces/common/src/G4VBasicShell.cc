// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VBasicShell.cc,v 1.2 1999-04-13 01:26:32 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VBasicShell.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4UImanager.hh"
#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

G4VBasicShell::G4VBasicShell()
{
  currentDirectory = "/";
}

G4VBasicShell::~G4VBasicShell() 
{;}

G4String G4VBasicShell::ModifyToFullPathCommand(const char* aCommandLine)
{
  G4String rawCommandLine = aCommandLine;
  if(rawCommandLine.isNull()||rawCommandLine(0)=='\0') return rawCommandLine;
  G4String commandLine = rawCommandLine.strip(G4String::both);
  G4String commandString;
  G4String parameterString;
  int i = commandLine.index(" ");
  if( i != RW_NPOS )
  {
    commandString = commandLine(0,i);
    parameterString = " ";
    parameterString += commandLine(i+1,commandLine.length()-(i+1));
  }
  else
  { commandString = commandLine; }

  G4String fullPathCommandLine
    = ModifyPath( commandString )+parameterString;
  return fullPathCommandLine;
}

G4String G4VBasicShell::GetCurrentWorkingDirectory()
{
  return currentDirectory;
}

G4bool G4VBasicShell::ChangeDirectory(const char* newDir)
{
  G4String aNewPrefix = newDir;
  G4String newPrefix = aNewPrefix.strip(G4String::both);
  G4String newDirectory = ModifyPath( newPrefix );
  if( newDirectory( newDirectory.length() - 1 ) != '/' )
  { newDirectory += "/"; }
  if( FindDirectory( newDirectory ) == NULL )
  { return false; }
  currentDirectory = newDirectory;
  return true;
}

G4UIcommandTree* G4VBasicShell::FindDirectory(const char* dirName)
{
  G4String aDirName = dirName;
  G4String theDir = aDirName.strip(G4String::both);
  G4String targetDir = ModifyPath( theDir );
  if( targetDir( targetDir.length()-1 ) != '/' )
  { targetDir += "/"; }
  G4UIcommandTree* comTree = G4UImanager::GetUIpointer()->GetTree();
  if( targetDir == "/" )
  { return comTree; }
  int idx = 1;
  while( idx < targetDir.length()-1 )
  {
    int i = targetDir.index("/",idx);
    comTree = comTree->GetTree(targetDir(0,i+1));
    if( comTree == NULL ) 
    { return NULL; }
    idx = i+1;
  }
  return comTree;
}

G4UIcommand* G4VBasicShell::FindCommand(const char* commandName)
{
  G4String rawCommandLine = commandName;
  G4String commandLine = rawCommandLine.strip(G4String::both);
  G4String commandString;
  int i = commandLine.index(" ");
  if( i != RW_NPOS )
  { commandString = commandLine(0,i); }
  else
  { commandString = commandLine; }

  G4String targetCom = ModifyPath(commandString);
  return G4UImanager::GetUIpointer()->GetTree()->FindPath(targetCom);
}

G4String G4VBasicShell::ModifyPath(G4String tempPath)
{
  G4String newPath = currentDirectory;
  if( tempPath(0) == '/' )   // full path is given
  { newPath = tempPath; }
  else if( tempPath(0) != '.' ) // add current prefix
  { newPath += tempPath; }
  else if( tempPath(0,2) == "./" ) // add current prefix
  { newPath += tempPath(2,tempPath.length()-2); }
  else                       // swim up with ".."
  {
    while( 1 )
    {
      if( tempPath(0,2) == ".." )
      {
        if( newPath != "/" )
        { 
	  G4String tmpString = newPath(0,newPath.length()-1);
          newPath = newPath(0,tmpString.last('/')+1); 
        }
        if( tempPath == ".." || tempPath == "../" )
        { break; }
        tempPath = tempPath(3,tempPath.length()-3);
      }
      else
      {
        newPath += tempPath;
        break;
      }
    }
  }
  return newPath;
}
////////////////////////////////////////////
// Method involving an interactive G4cout //
////////////////////////////////////////////
/***************************************************************************/
void G4VBasicShell::ExecuteCommand (
 G4String aCommand
)
/***************************************************************************/
// Should be put in G4VBasicShell.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(aCommand.length()<2) return;
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  int commandStatus = UI->ApplyCommand(aCommand);
  switch(commandStatus) {
  case fCommandSucceeded:
    break;
  case fCommandNotFound:
    G4cerr << "command not found" << endl;
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
/***************************************************************************/
void G4VBasicShell::ApplyShellCommand (
 G4String a_string
,G4bool& exitSession
,G4bool& exitPause
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;

  G4String command = a_string.strip(G4String::leading);
  if( command(0) == '#' ) { 

    G4cout << command << endl; 

  } else if( command == "ls" || command(0,3) == "ls " ) {

    ListDirectory( command );

  } else if( command == "pwd" ) { 

    G4cout << "Current Working Directory : " 
       << GetCurrentWorkingDirectory() << endl; 

  } else if( command(0,2) == "cd" ) { 

    ChangeDirectoryCommand ( command );

  } else if( command(0,4) == "help" ) { 

    TerminalHelp( command ); 

  } else if( command(0) == '?' ) { 

    ShowCurrent( command );

  } else if( command(0,4) == "hist" ) {

    G4int nh = UI->GetNumberOfHistory();
    for(int i=0;i<nh;i++) { 
      G4cout << i << ": " << UI->GetPreviousCommand(i) << endl; 
    }

  } else if( command(0) == '!' ) {

    G4String ss = command(1,command.length()-1);
    G4int vl;
    const char* tt = ss;
    istrstream is((char*)tt);
    is >> vl;
    G4int nh = UI->GetNumberOfHistory();
    if(vl>=0 && vl<nh) { 
      G4String prev = UI->GetPreviousCommand(vl); 
      G4cout << prev << endl;
      ExecuteCommand (ModifyToFullPathCommand(prev));
    } else { 
      G4cerr << "history " << vl << " is not found." << endl; 
    }

  } else if( command(0,4) == "exit" ) { 

    if( exitPause == false) { //In a secondary loop.
      G4cout << "You are now processing RUN." << endl;
      G4cout << "Please abort it using \"/run/abort\" command first" << endl;
      G4cout << " and use \"continue\" command until the application" << endl;
      G4cout << " becomes to Idle." << endl;
    } else {
      exitSession = true;
    }

  } else if( command(0,4) == "cont" ) { 

    exitPause = true;

  } else {

    ExecuteCommand (ModifyToFullPathCommand(a_string));

  }
}
void G4VBasicShell::ShowCurrent(G4String newCommand)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
  G4String comString = newCommand(1,newCommand.length()-1);
  G4String theCommand = ModifyToFullPathCommand(comString);
  G4String curV = UI->GetCurrentValues(theCommand);
  if( ! curV.isNull() ) { 
    G4cout << "Current value(s) of the parameter(s) : " << curV << endl; 
  }
}
void G4VBasicShell::ChangeDirectoryCommand(G4String newCommand)
{
  G4String prefix;
  if( newCommand.length() <= 3 ) { 
    prefix = "/"; 
  } else {
    G4String aNewPrefix = newCommand(3,newCommand.length()-3);
    prefix = aNewPrefix.strip(G4String::both);
  }
  if(!ChangeDirectory(prefix)) { 
    G4cout << "directory <" << prefix << "> not found." << endl; 
  }
}
void G4VBasicShell::ListDirectory(G4String newCommand)
{
  G4String targetDir;
  if( newCommand.length() <= 3 ) { 
    targetDir = GetCurrentWorkingDirectory();
  } else {
    G4String newPrefix = newCommand(3,newCommand.length()-3);
    targetDir = newPrefix.strip(G4String::both);
  }
  G4UIcommandTree* commandTree = FindDirectory( targetDir );
  if( commandTree == NULL ) { 
    G4cout << "Directory <" << targetDir << "> is not found." << endl; 
  } else { 
    commandTree->ListCurrent(); 
  }
}
void G4VBasicShell::TerminalHelp(G4String newCommand)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;
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
   //G4cout << endl << "Type the number ( 0:end, -n:n level back ) : "<<flush;
    G4cout << endl << "Type the number ( 0:end, -n:n level back ) : "<<endl;
    G4int i;
    if(!GetHelpChoice(i)){
      G4cout << endl << "Not a number, once more" << endl; 
      continue;
    } else if( i < 0 ){
      iFloor += i;
      if( iFloor < 0 ) iFloor = 0;
      floor[iFloor]->ListCurrentWithNum(); 
      continue;
    } else if(i == 0) { 
      break;
    } else if( i > 0 ) {
      int n_tree = floor[iFloor]->GetTreeEntry();
      if( i > n_tree )
      { 
        if( i <= n_tree + floor[iFloor]->GetCommandEntry() )
        { 
          floor[iFloor]->GetCommand(i-n_tree)->List(); 
        }
      }
      else
      {
        floor[iFloor+1] = floor[iFloor]->GetTree(i);
        iFloor++;
        floor[iFloor]->ListCurrentWithNum();
      }
    }
  }
  G4cout << "Exit from HELP." << endl << endl;
  //G4cout << endl;
  ExitHelp();
}


