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
// $Id$
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
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout <<
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    "\n!!!!! Xaw is deprecated and will be removed in the next major release."
    "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    << G4endl;
  }

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  G4Xt*     interactorManager = G4Xt::getInstance (argc,argv,(char*)"Xaw");
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
const G4String& aPrompt
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
const G4String& a_state
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
) const
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
void G4UIXaw::Callback (
 Widget
,XtPointer a_tag
,XtPointer
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

  //a_widget     = NULL;  Not used (1st argument).  Comment out to avoid compiler warnings. (JA)
  //a_data       = NULL;  Not used (3rd argument).  Comment out to avoid compiler warnings. (JA)
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
