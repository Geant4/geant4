// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXaw.cc,v 1.1 1999-01-07 16:09:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

//#define DEBUG

#ifdef G4UI_BUILD_XAW_SESSION

#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <Xaw/Dialog.h>
#include <Xaw/Command.h>

#include "G4UIXaw.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommandStatus.hh"
#include "G4Xt.hh"


static void  Callback                        (Widget,XtPointer,XtPointer);

static G4bool exitSession = true;
static G4bool exitPause   = true;
/***************************************************************************/
G4UIXaw::G4UIXaw (
 int argc
,char** argv
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  G4Xt*     interactorManager = G4Xt::getInstance (argc,argv,"Xaw");
  Widget    top = (Widget)interactorManager->GetMainInteractor();

  shell     = XtAppCreateShell      ("G4UIXaw","G4UIXaw",topLevelShellWidgetClass,XtDisplay(top),NULL,0); 

  Arg       args[2];
  XtSetArg  (args[0],XtNlabel,"G4 command");
  XtSetArg  (args[1],XtNvalue,"");             // Needed to have a text Area.
  dialog    = XtCreateManagedWidget ("dialog",dialogWidgetClass,shell,args,2);

  XawDialogAddButton (dialog,"Ok",Callback,(XtPointer)this);

  XtRealizeWidget (shell);
}
/***************************************************************************/
G4UIXaw::~G4UIXaw (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{ 
  XtDestroyWidget (shell);
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(NULL);
}
/***************************************************************************/
G4UIsession* G4UIXaw::SessionStart (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xt*        interactorManager = G4Xt::getInstance ();
  Prompt       ("session");
  exitSession  = false;
  interactorManager->DisableSecondaryLoop ();
  void*        event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitSession==true) break;
  }
  interactorManager->EnableSecondaryLoop ();
  return       this;
}
/***************************************************************************/
void G4UIXaw::Prompt (
 G4String aPrompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  Arg          args[1];
  XtSetArg     (args[0],XtNlabel,aPrompt.data());
  XtSetValues  (dialog,args,1);
}
/***************************************************************************/
void G4UIXaw::SessionTerminate (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIXaw::PauseSessionStart (
 G4String a_state
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_state=="G4_pause> ")
    { 
      SecondaryLoop ("Pause, type continue to exit this state");
    }

  if(a_state=="EndOfEvent")
    {
      // Picking with feed back in event data Done here !!!
      SecondaryLoop ("End of event, type continue to exit this state");
    }
}
/***************************************************************************/
void G4UIXaw::SecondaryLoop (
 G4String a_prompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xt*         interactorManager = G4Xt::getInstance ();
  Prompt        (a_prompt);
  exitPause     = false;
  void*         event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitPause==true) break;
  }
  Prompt       ("session");
}
/***************************************************************************/
void G4UIXaw::ApplyShellCommand (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI==NULL) return;

  G4String     command = a_string.strip(G4String::leading);
  if( command(0) == '#' ) { 

    G4cout << command << endl; 

  } else if( command(0,2) == "ls" ) { 

    ListDirectory( command );

  } else if( command == "pwd" ) { 

    G4cout << "Current Working Directory : " 
       << GetCurrentWorkingDirectory() << endl; 

  } else if( command(0,2) == "cd" ) { 

    ChangeDirectoryCommand ( command );

  } else if( command(0,4) == "help" ) { 

    //TerminalHelp( command ); 
    G4cout << "Not implemented." << endl; 

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
/***************************************************************************/
void G4UIXaw::ExecuteCommand (
 G4String aCommand
)
/***************************************************************************/
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
void G4UIXaw::ShowCurrent ( 
 G4String newCommand 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
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
/***************************************************************************/
void G4UIXaw::ChangeDirectoryCommand ( 
 G4String newCommand 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
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
/***************************************************************************/
void G4UIXaw::ListDirectory( 
 G4String newCommand 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4String targetDir;
  if( newCommand.length() <= 3 ) { 
    targetDir = "./"; 
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
/***************************************************************************/
Widget G4UIXaw::GetDialog (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return dialog;
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
static void Callback (
 Widget a_widget
,XtPointer a_tag
,XtPointer a_data
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UIXaw*     This = (G4UIXaw*)a_tag;
  Widget       dialog = This->GetDialog();
  char*        value = XawDialogGetValueString(dialog);
  if(value==NULL) return;
  G4String     command (value);

  This->ApplyShellCommand (command);

  Arg          args[1];
  XtSetArg     (args[0],XtNvalue,"");
  XtSetValues  (dialog,args,1);

  a_widget     = NULL;
  a_data       = NULL;
}

#endif








