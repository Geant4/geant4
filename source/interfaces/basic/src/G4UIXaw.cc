// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIXaw.cc,v 1.2 1999-04-13 01:26:27 yhajime Exp $
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

static G4bool ConvertStringToInt(const char*,int&);

static G4bool exitSession = true;
static G4bool exitPause = true;
static G4bool exitHelp = true;
/***************************************************************************/
G4UIXaw::G4UIXaw (
 int argc
,char** argv
)
:fHelp(false)
,fHelpChoice(0)
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
  if(a_state=="G4_pause> ") { 
    SecondaryLoop ("Pause, type continue to exit this state");
  }

  if(a_state=="EndOfEvent") {
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
Widget G4UIXaw::GetDialog (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return dialog;
}
/***************************************************************************/
G4bool G4UIXaw::GetHelpChoice(
 G4int& aInt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  fHelp = true;
  //
  G4Xt* interactorManager = G4Xt::getInstance ();
  Prompt("Help");
  exitHelp = false;
  void* event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    interactorManager->DispatchEvent(event);
    if(exitHelp==true) break;
  }
  Prompt("session");
  //
  if(fHelp==false) return false;
  aInt = fHelpChoice;
  fHelp = false;
  return true;
}
/***************************************************************************/
void G4UIXaw::ExitHelp(
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void G4UIXaw::Callback (
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

  if(This->fHelp==true) {
    exitHelp = true;
    This->fHelp = ConvertStringToInt(command.data(),This->fHelpChoice);
  } else {
    This->ApplyShellCommand (command,exitSession,exitPause);
  }

  Arg          args[1];
  XtSetArg     (args[0],XtNvalue,"");
  XtSetValues  (dialog,args,1);

  a_widget     = NULL;
  a_data       = NULL;
}
//////////////////////////////////////////////////////////////////////////////
G4bool ConvertStringToInt(
 const char* aString
,int& aInt
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  aInt = 0;
  if(aString==NULL) return false;
  char* s;
  long value = strtol(aString,&s,10);
  if(s==aString) return false;
  aInt = value;
  return true;
}

#endif
