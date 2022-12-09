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
//
// G.Barrand

//#define DEBUG

#include "G4Types.hh"

#include <string.h>

#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <X11/keysym.h>

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

#include <stdlib.h>

static void XmTextAppendString (Widget,char*);

static void clearButtonCallback (Widget,XtPointer,XtPointer);

static char* XmConvertCompoundStringToString (XmString,G4int);
static G4bool ConvertStringToInt(const char*,G4int&);
static void ExecuteChangeSizeFunction(Widget);

static G4bool exitSession = true;
static G4bool exitPause = true;
static G4bool exitHelp = true;
/***************************************************************************/
G4UIXm::G4UIXm (
 G4int argc
,char** argv
)
:shell(NULL)
,command(NULL)
,menuBar(NULL)
,text(NULL)
,fHelp(false)
,fHelpChoice(0)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  G4Xt* interactorManager = G4Xt::getInstance (argc,argv,(char*)"Xm");

  Widget top = (Widget)interactorManager->GetMainInteractor();

  if(std::getenv("XENVIRONMENT")==NULL) {
    XrmDatabase database = XrmGetDatabase(XtDisplay(top));
    if(database!=NULL) {
      XrmPutLineResource(&database,"*topShadowColor:white");
      XrmPutLineResource(&database,"*bottomShadowColor:black");
      XrmPutLineResource(&database,"*foreground:black");
      XrmPutLineResource(&database,"*background:lightgrey");
      XrmPutLineResource(&database,"*borderColor:lightgrey");
      XrmPutLineResource(&database,"*fontList:-*-helvetica-bold-r-*-*-*-120-*-*-*-*-iso8859-1");
      XrmPutLineResource(&database,"*text.background:white");
      XrmPutLineResource(&database,"*text.fontList:*courier*-r-*--14-*");
      XrmPutLineResource(&database,"*text.maxLength:8000");
    }
  }

  Arg args[9];
  XtSetArg(args[0],XmNkeyboardFocusPolicy,XmPOINTER); // For completion.
  shell = XtAppCreateShell ("G4UIXm","G4UIXm",
			    topLevelShellWidgetClass,XtDisplay(top),
			    args,1);
  form = XmCreateForm (shell,(char*)menu_str[0].c_str(),NULL,0);
  XtManageChild (form);

  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_FORM);
  XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
  menuBar = XmCreateMenuBar (form,(char*)menu_str[1].c_str(),args,3);

  XtSetArg(args[0],XmNtopAttachment      ,XmATTACH_NONE);
  XtSetArg(args[1],XmNleftAttachment     ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment    ,XmATTACH_FORM);
  XtSetArg(args[3],XmNbottomAttachment   ,XmATTACH_FORM);
  command = XmCreateCommand (form,(char*)menu_str[2].c_str(),args,4);
  XtManageChild (command);

  XtSetArg(args[0],XmNtopAttachment   ,XmATTACH_NONE);
  XtSetArg(args[1],XmNleftAttachment  ,XmATTACH_FORM);
  XtSetArg(args[2],XmNrightAttachment ,XmATTACH_FORM);
  XtSetArg(args[3],XmNbottomAttachment,XmATTACH_WIDGET);
  XtSetArg(args[4],XmNbottomWidget    ,command);
  XmString cps = XmStringLtoRCreate((char*)menu_str[3].c_str(),(char*)XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[5],XmNlabelString,cps);
  Widget clearButton = XmCreatePushButton(form,(char*)menu_str[4].c_str(),args,6);
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
  text = XmCreateScrolledText (form,(char*)menu_str[5].c_str(),args,9);
  XtManageChild (text);

  XtAddCallback(clearButton,XmNactivateCallback,
		clearButtonCallback,(XtPointer)text);
  XtAddCallback(command,XmNcommandEnteredCallback,
		CommandEnteredCallback,(XtPointer)this);

  Widget commandText = XmCommandGetChild(command,XmDIALOG_COMMAND_TEXT);
  XtAddEventHandler(commandText,KeyPressMask,False,keyHandler,(XtPointer)this);

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
  XmString cps = XmStringLtoRCreate(str,(char*)XmSTRING_DEFAULT_CHARSET);
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
G4int G4UIXm::ReceiveG4cout (
 const G4String& a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XmTextAppendString(text,(char*)a_string.data());
  return 0;
}
/***************************************************************************/
G4int G4UIXm::ReceiveG4cerr (
 const G4String& a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  XmTextAppendString(text,(char*)a_string.data());
  return 0;
}
/***************************************************************************/
G4bool G4UIXm::GetHelpChoice(
 G4int& aInt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  fHelp = true;
  // SecondaryLoop :
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
void G4UIXm::ExitHelp(
) const
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
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
  XmString cps = XmStringLtoRCreate((char*)a_label,(char*)XmSTRING_DEFAULT_CHARSET);
  XtSetArg (args[0],XmNlabelString,cps);
  XtSetArg (args[1],XmNsubMenuId,widget);
  widget = XmCreateCascadeButton (menuBar,(char*)a_name,args,2);
  XmStringFree (cps);
  XtManageChild (widget);
  ExecuteChangeSizeFunction(form);
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
  XtAddCallback (widget,XmNactivateCallback,ButtonCallback,(XtPointer)this);
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
void G4UIXm::CommandEnteredCallback (
 Widget
,XtPointer a_tag
,XtPointer a_data
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UIXm*  This = (G4UIXm*)a_tag;

  XmString cps  = ((XmCommandCallbackStruct*)a_data)->value;
  char*    ss = XmConvertCompoundStringToString(cps,0);
  G4String scommand (ss);
  XtFree   (ss);

  if(This->fHelp==true) {
    exitHelp = true;
    This->fHelp = ConvertStringToInt(scommand.data(),This->fHelpChoice);
  } else {
    This->ApplyShellCommand (scommand,exitSession,exitPause);
  }

  a_tag    = NULL;
}
/***************************************************************************/
void G4UIXm::keyHandler (
 Widget a_widget
,XtPointer a_tag
,XEvent* a_event
,Boolean*
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  KeySym keySym;
  XLookupString(&(a_event->xkey),NULL,0,&keySym,NULL);
  if(keySym!=XK_Tab) return;
  G4UIXm* This = (G4UIXm*)a_tag;
  char* s = XmTextGetString(a_widget);
  G4String ss = This->Complete(s);
  XmTextSetString(a_widget,(char*)ss.data());
  XtFree(s);
  XmTextSetInsertionPosition(a_widget,XmTextGetLastPosition(a_widget));
}
/***************************************************************************/
void clearButtonCallback (
 Widget
,XtPointer a_tag
,XtPointer
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  static const G4String empty_str = "";
  XmTextSetString((Widget)a_tag,(char*)empty_str.c_str());
}
/***************************************************************************/
void G4UIXm::ButtonCallback (
 Widget a_widget
,XtPointer a_tag
,XtPointer
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UIXm* This = (G4UIXm*)a_tag;
  if(This->fHelp==true) return; // Disabled when in help.
  G4String ss = This->GetCommand (a_widget);
  This->ApplyShellCommand(ss,exitSession,exitPause);
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
char* XmConvertCompoundStringToString (
 XmString a_cps 
,G4int a_number 
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_cps==NULL) return NULL;
  char* ss = NULL;
  XmStringContext context;
  XmStringInitContext(&context,a_cps);
  G4int icount = 0;
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
//////////////////////////////////////////////////////////////////////////////
G4bool ConvertStringToInt(
 const char* aString
,G4int& aInt
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  aInt = 0;
  if(aString==NULL) return false;
  char* s;
  G4long value = strtol(aString,&s,10);
  if(s==aString) return false;
  aInt = (G4int)value;
  return true;
}
#include <X11/IntrinsicP.h>
//////////////////////////////////////////////////////////////////////////////
void ExecuteChangeSizeFunction (
 Widget aWidget
)
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(aWidget==NULL) return;
  if(aWidget->core.widget_class->core_class.resize==NULL) return;
  (aWidget->core.widget_class->core_class.resize)(aWidget);
}
