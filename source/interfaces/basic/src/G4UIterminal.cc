// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIterminal.cc,v 1.5 1999-11-15 10:39:44 gunter Exp $
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
       if ( aCommand.index("@@") != G4std::string::npos) {
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
    newCommand.readLine( cin, FALSE );
    if (!cin.good()) { cin.clear(); newCommand = nullString; iExit=false; break; }
    newCommand = newCommand.strip(G4String::leading);
    if( newCommand.length() < 1 ) { break; }

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
    if( nC.length() < 1 ){ break; }
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

G4bool G4UIterminal::GetHelpChoice(G4int& aInt){
  cin >> aInt;
  if(!cin.good()){
    cin.clear();
    cin.ignore(30,'\n');
    return false;
  } 
  return true;
}
void G4UIterminal::ExitHelp(){
  // cin.flush();
  char temp[100];
  cin.getline( temp, 100 );
}




