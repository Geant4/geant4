// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXm.cc,v 1.1 1999-01-07 16:09:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G.Barrand

//#define DEBUG

#ifdef G4UI_BUILD_XM_SESSION

#ifdef WIN32
#include <Strstrea.h>
#else
#include <strstream.h>
#endif

#include <X11/Intrinsic.h>
#include <X11/Shell.h>

#include <Xm/Xm.h>
#include <Xm/Command.h>
#include <Xm/RowColumn.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/CascadeB.h>
#include <Xm/Text.h>

#include "G4UIXm.hh"
#include "G4UImanager.hh"
#include "G4StateManager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommandStatus.hh"

#include "G4Xt.hh"

static void XmTextAppendString (Widget,char*);

static void commandCallback (Widget,XtPointer,XtPointer);
static void buttonCallback (Widget,XtPointer,XtPointer);
static void clearButtonCallback (Widget,XtPointer,XtPointer);
static char* XmConvertCompoundStringToString (XmString,int);

static G4bool exitSession = true;
static G4bool exitPause   = true;
static unsigned CommandsHashFun (const Widget& interactor) {
  // Is it correct to return a pointer ?
  return (unsigned)(unsigned long)interactor;
}
/***************************************************************************/
G4UIXm::G4UIXm (
 int argc
,char** argv
)
:shell(NULL)
,command(NULL)
,menuBar(NULL)
,text(NULL)
,commands(CommandsHashFun)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  G4Xt* interactorManager = G4Xt::getInstance (argc,argv,"Xm");

  Widget top = (Widget)interactorManager->GetMainInteractor();

  shell = XtAppCreateShell ("G4UIXm","G4UIXm",
			    topLevelShellWidgetClass,XtDisplay(top),
			    NULL,0); 
  Widget form = XmCreateForm (shell,"form",NULL,0);
  XtManageChild (form);

  Arg args[9];
  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
  XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
  menuBar = XmCreateMenuBar (form,"menuBar",args,3);

  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_NONE);
  XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
  XtSetArg(args[3],XmNbottomAttachment,XmATTACH_FORM);
  command = XmCreateCommand (form,"command",args,4);
  XtManageChild (command);

  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_NONE);
  XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
  XtSetArg(args[3],XmNbottomAttachment,XmATTACH_WIDGET);
  XtSetArg(args[4],XmNbottomWidget    ,command);
  XmString cps = XmStringLtoRCreate("Clear",XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[5],XmNlabelString,cps);
  Widget clearButton = XmCreatePushButton(form,"clearButton",args,6);
  XmStringFree (cps);
  XtManageChild (clearButton);

  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_WIDGET);
  XtSetArg(args[1],XmNtopWidget       ,menuBar);
  XtSetArg(args[2],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[3],XmNrightAttachment ,XmATTACH_FORM);
  XtSetArg(args[4],XmNbottomAttachment,XmATTACH_WIDGET);
  XtSetArg(args[5],XmNbottomWidget    ,clearButton);
  XtSetArg(args[6],XmNeditMode        ,XmMULTI_LINE_EDIT);
  XtSetArg(args[7],XmNrows            ,12);
  XtSetArg(args[8],XmNcolumns         ,80);
  text = XmCreateScrolledText (form,"text",args,9);
  XtManageChild (text);

  XtAddCallback(clearButton,XmNactivateCallback,
		clearButtonCallback,(XtPointer)text);
  XtAddCallback(command,XmNcommandEnteredCallback,
		commandCallback,(XtPointer)this);

  XtRealizeWidget(shell);
  XtMapWidget(shell);

  if(UI!=NULL) UI->SetCoutDestination(this);
}
/***************************************************************************/
G4UIXm::~G4UIXm(
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{ 
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) {
    UI->SetSession(NULL);
    UI->SetCoutDestination(NULL);
  }
  XtDestroyWidget(shell);
}
/***************************************************************************/
G4UIsession* G4UIXm::SessionStart (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xt* interactorManager = G4Xt::getInstance ();
  Prompt("session");
  exitSession = false;
  interactorManager->DisableSecondaryLoop ();
  void* event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitSession==true) break;
  }
  interactorManager->EnableSecondaryLoop ();
  return this;
}
/***************************************************************************/
void G4UIXm::Prompt (
 G4String aPrompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  Arg args[1];
  char* str = (char*)XtNewString(aPrompt.data());
  XmString cps = XmStringLtoRCreate(str,XmSTRING_DEFAULT_CHARSET);
  XtFree(str);
  XtSetArg(args[0],XmNpromptString,cps);
  XtSetValues(command,args,1);
  XmStringFree(cps);
}
/***************************************************************************/
void G4UIXm::SessionTerminate (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIXm::PauseSessionStart (
 G4String a_state
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_state=="G4_pause> ") { 
    SecondaryLoop ("Pause, type continue to exit this state");
  }

  if(a_state=="EndOfEvent") {
    // Picking with feed back in event data Done here !!!
    SecondaryLoop ("End of event, type continue to exit this state");
  }
}
/***************************************************************************/
void G4UIXm::SecondaryLoop (
 G4String a_prompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xt* interactorManager = G4Xt::getInstance ();
  Prompt(a_prompt);
  exitPause = false;
  void* event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitPause==true) break;
  }
  Prompt("session");
}
/***************************************************************************/
void G4UIXm::ApplyShellCommand (
 G4String a_string
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
void G4UIXm::ExecuteCommand (
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
G4int G4UIXm::ReceiveG4cout (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XmTextAppendString(text,(char*)a_string.data());
  return 0;
}
/***************************************************************************/
G4int G4UIXm::ReceiveG4cerr (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XmTextAppendString(text,(char*)a_string.data());
  return 0;
}
/***************************************************************************/
void G4UIXm::ShowCurrent ( 
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
void G4UIXm::ChangeDirectoryCommand ( 
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
void G4UIXm::ListDirectory( 
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
void G4UIXm::AddMenu (
 const char* a_name
,const char* a_label
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(menuBar==NULL) return;
  if(a_name==NULL) return;
  if(a_label==NULL) return;
  XtManageChild (menuBar);
  // Pulldown menu :
  Widget widget;
  widget = XmCreatePulldownMenu (menuBar,(char*)a_name,NULL,0);
  AddInteractor (a_name,(G4Interactor)widget);
  // Cascade button :
  Arg args[2];
  XmString cps = XmStringLtoRCreate((char*)a_label,XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[0],XmNlabelString,cps);
  XtSetArg (args[1],XmNsubMenuId,widget);
  widget = XmCreateCascadeButton (menuBar,(char*)a_name,args,2);
  XmStringFree (cps);
  XtManageChild (widget);
}
/***************************************************************************/
void G4UIXm::AddButton (
 const char* a_menu
,const char* a_label
,const char* a_command
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_menu==NULL) return;
  if(a_label==NULL) return;
  if(a_command==NULL) return;
  Widget parent = (Widget)GetInteractor(a_menu);
  if(parent==NULL) return;
  Widget widget = XmCreatePushButton(parent,(char*)a_label,NULL,0);
  XtManageChild (widget);
  XtAddCallback (widget,XmNactivateCallback,buttonCallback,(XtPointer)this);
  commands[widget] = a_command;
}
/***************************************************************************/
G4String G4UIXm::GetCommand (
 Widget a_widget
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return commands[a_widget];
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void commandCallback (
 Widget    a_widget
,XtPointer a_tag
,XtPointer a_data
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UIXm*  This = (G4UIXm*)a_tag;

  XmString cps  = ((XmCommandCallbackStruct*)a_data)->value;
  char*    ss = XmConvertCompoundStringToString(cps,0);
  G4String command (ss);
  XtFree   (ss);

  This->ApplyShellCommand (command);

  a_widget = NULL;
  a_tag    = NULL;
}
/***************************************************************************/
void clearButtonCallback (
 Widget a_widget
,XtPointer a_tag
,XtPointer
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XmTextSetString((Widget)a_tag,"");
}
/***************************************************************************/
void buttonCallback (
 Widget a_widget
,XtPointer a_tag
,XtPointer
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UIXm* This = (G4UIXm*)a_tag;
  G4String ss = This->GetCommand (a_widget);
  //printf ("debug : execute:\n%s\n",ss.data());
  This->ApplyShellCommand (ss);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
char* XmConvertCompoundStringToString (
 XmString a_cps 
,int a_number 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_cps==NULL) return NULL;
  char* ss = NULL;
  XmStringContext context;
  XmStringInitContext(&context,a_cps);
  int icount = 0;
  Boolean Done = False;
  while(Done==False) {  
    char* text = NULL;
    XmStringCharSet charset = NULL;
    XmStringDirection direct;
    Boolean sep;
    if(XmStringGetNextSegment(context,&text,&charset,&direct,&sep)==True) {
      XtFree(charset);
      if(sep==True) Done = True;
      if(icount==a_number) { 
	ss = text;
	break;
      }
      icount++;
      XtFree(text);
    }
    else
      Done = True;
  }
  XmStringFreeContext(context);
  return ss;
}
/***************************************************************************/
void XmTextAppendString (
 Widget This
,char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(This==NULL) return;
  if(!XtIsSubclass(This,xmTextWidgetClass)) return;
  if(a_string==NULL) return;
  XmTextPosition  lastpos = XmTextGetLastPosition(This);
  XmTextReplace(This,lastpos,lastpos,a_string);
  XmTextSetInsertionPosition(This,XmTextGetLastPosition(This));
}

#endif
