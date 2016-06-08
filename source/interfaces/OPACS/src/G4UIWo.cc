// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIWo.cc,v 1.4 2000/08/11 09:56:50 barrand Exp $
// GEANT4 tag $Name: geant4-03-01 $
//
//G4UIWo.cc   -- copied from G4UIterminal.cc

//#define DEBUG

#ifdef G4UI_BUILD_WO_SESSION

#include <stdlib.h>

#ifdef HAS_XM
#include <Xm/MainW.h> //For G4OIX
#include <Xm/Label.h>
#endif
#ifdef HAS_XAW
#include <X11/Intrinsic.h>
#include <Xaw/Label.h>
#endif

// Below are the "wanted" packages header files.
// begin Want_h
#include <WoWo.h>
#include <WoXtw.h>
#include <WoXm.h>
#include <WoXaw.h>
#include <WoXo.h>
//end Want_h

/*Co*/
#include <CPrinter.h>
#include <CString.h>
#include <CFile.h>
#include <CText.h>
#include <OShell.h>
/*Xx*/
#include <XWidget.h>
/*Xo*/
#include <XoCamera.h>
/*Wo*/
#include <OEvent.h>
#include <OWidget.h>
#include <OInterface.h>
/*G4*/
#include "G4UIWo.hh"
#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4WoMessenger.hh"
#include "G4StateManager.hh"
#include "G4Xt.hh"
#include "G4oCommon.hh"

static Widget CreateG4View    (String,Widget,String,ArgList,Cardinal);
static Widget CreateG4XoView  (Widget,String,ArgList,Cardinal);
static Widget CreateG4OIXView (Widget,String,ArgList,Cardinal);

static int Inited = 0;
/***************************************************************************/
G4UIWo::G4UIWo (
 int argc
,char** argv
)
:woMessenger(NULL)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) UI->SetSession(this);

  if(Inited==1) return; // Below code must be executed once.

// Below are the "wanted" packages include c files.
// begin Want_c
#include <WoWo.ic>
#include <WoXtw.ic>
#include <WoXm.ic>
#include <WoXaw.ic>
#include <WoXo.ic>
// end Want_c

  //Declare G4 interpreter to Wo.
  OInterpreterCreate  ("G4",(OInterpreterDoFunction)ExecuteScript,NULL,NULL,NULL);
  OShellAddNewCommand (WoGetShellInterpreter(),"G4/G4",Execute_G4);
  SetTypes            ();

  //Declare osh to G4.
  woMessenger = new G4WoMessenger (WoGetShellInterpreter());

  //Declare G4/vis~/views to Wo.
#ifdef HAS_XO
  OClassDeclareCompoundWidget ("G4/G4XoView" ,(OClassCreateWidgetFunction)CreateG4XoView ,xoCameraWidgetClass);
#endif
#ifdef HAS_XM
  OClassDeclareCompoundWidget ("G4/G4OIXView",(OClassCreateWidgetFunction)CreateG4OIXView,xmMainWindowWidgetClass);
  OClassDeclareWidgetClass    ("G4/ThumbWheel"       ,NULL);
  OClassDeclareWidgetClass    ("G4/SoGLwMDrawingArea",NULL);
  OClassDeclareWidgetClass    ("G4/XoSoArea"         ,NULL);
#endif

  G4Xt*        interactorManager = G4Xt::getInstance (argc,argv,"Wo");
  Widget       top = (Widget)interactorManager->GetMainInteractor();
  XtAppContext context = XtWidgetToApplicationContext(top);

  WoSetXt      (context,top);
  interactorManager->AddDispatcher ((G4DispatchFunction)WoDispatchEvent);

  // Load file in SessionStart.
  // At this moment the G4VisManager is expected to be created.
  // It is needed if an odb file is expected to create some upper G4View.
  WoSetInterfaceFile ("nil"); 
  WoStartup          (argc,argv);

  if(UI!=NULL) UI->SetCoutDestination(this);

  Inited = 1;
}
/***************************************************************************/
G4UIWo::~G4UIWo (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{ 
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) {
    UI->SetSession(NULL);
    UI->SetCoutDestination(NULL);
  }
  if(woMessenger!=NULL) delete woMessenger;
}
/***************************************************************************/
G4UIsession* G4UIWo::SessionStart (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  static int Done = 0;
  G4Xt* interactorManager = G4Xt::getInstance ();
  if(Done==0) { //Load odb file.
    char* fileName = NULL;
    if(getenv("WOENVIRONMENT")!=NULL) {
      fileName = CStringDuplicate (getenv("WOENVIRONMENT"));
    } else {
      char**     argv;
      int        argc;
      argv       = interactorManager->GetArguments (&argc);
      if(argc>=1) {
	char*         aname = CFileGetName   (argv[0]);
	fileName      = CStringCreateF (strlen(aname)+4,"%s.odb",aname);
	CStringDelete (aname);
      }
    }
    WoSetInterfaceFile             (fileName);
    OInterfaceLoadFile             (fileName,NULL);
    OInterfaceSetExtentNotModified ();
    if(OWidgetHasMappedShell()==0) XWidgetMap (OWidgetGetFirstShell());
    CStringDelete                  (fileName);
    Done = 1;
  }
  int          value;
  //  Prompt       ("session");
  interactorManager->DisableSecondaryLoop ();
  void*         event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    if(OEventIsExit((XEvent*)event,&value)==1) 
      {
	WoDispatchEvent (event);
	break;
      }
    interactorManager->DispatchEvent(event);
  }
  interactorManager->EnableSecondaryLoop ();
  return        this;
}
/***************************************************************************/
void G4UIWo::Prompt (
 G4String aPrompt
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  CWarnF("%s\n",aPrompt.data());
}
/***************************************************************************/
void G4UIWo::SessionTerminate (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4UIWo::PauseSessionStart (
 G4String a_state
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "debug : G4UIWo::PauseSessionStart : begin : " << a_state << G4endl;
#endif

  if(a_state=="G4_pause> ")
    { 
      SecondaryLoop ("G4Pause.osh");
    }

  if(a_state=="BeginOfRun")
    {
      WoExecuteShellFileInSameContext ("G4RunBegin.osh");
    }
  else if(a_state=="EndOfRun")
    {
      WoExecuteShellFileInSameContext ("G4RunEnd.osh");
    }
  else if(a_state=="BeginOfEvent")
    {
      WoExecuteShellFileInSameContext ("G4EventBegin.osh");
    }
  else if(a_state=="EndOfEvent")
    {
      // Picking with feed back in event data Done here !!!
      SecondaryLoop ("G4EventEnd.osh");
    }

#ifdef DEBUG
  G4cout << "debug : G4UIWo::PauseSessionStart : end." << G4endl;
#endif
}
/***************************************************************************/
void G4UIWo::SecondaryLoop (
 G4String a_file
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xt*        interactorManager = G4Xt::getInstance ();
  int          value;
  //  Prompt       (a_prompt);
  WoExecuteShellFileInSameContext ((char*)a_file.data());
  void*         event;
  while((event = interactorManager->GetEvent())!=NULL) { 
    if(OEventIsExit((XEvent*)event,&value)==1) {
      WoDispatchEvent (event);
      break;
    }
    interactorManager->DispatchEvent(event);
  }
  //  Prompt       ("session");
  WoExecuteShellFileInSameContext ((char*)"G4Session.osh");
}
/***************************************************************************/
G4int G4UIWo::ReceiveG4cout (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  CWarnF ("%s",a_string.data());
  return 0;
}
/***************************************************************************/
G4int G4UIWo::ReceiveG4cerr (
 G4String a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  CWarnF ("%s",a_string.data());
  return 0;
}
/***************************************************************************/
Widget CreateG4XoView (
 Widget   a_parent
,String   a_name
,ArgList  a_list
,Cardinal a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  return CreateG4View("Xo",a_parent,a_name,a_list,a_number);
}
/***************************************************************************/
Widget CreateG4OIXView (
 Widget   a_parent
,String   a_name
,ArgList  a_list
,Cardinal a_number
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  return CreateG4View("OIX",a_parent,a_name,a_list,a_number);
}
/***************************************************************************/
Widget CreateG4View (
 String   a_gs
,Widget   a_parent
,String   a_name
,ArgList  /*a_list*/
,Cardinal /*a_number*/
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  Widget widget = NULL;
/*.........................................................................*/
  if(a_gs==NULL)     return NULL;
  if(a_parent==NULL) return NULL;
  if(a_name==NULL)   return NULL;

  G4UImanager* UI = G4UImanager::GetUIpointer();
  if(UI!=NULL) {

    G4Xt*  interactorManager     = G4Xt::getInstance ();

    interactorManager->SetParentInteractor (a_parent);
    interactorManager->SetCreationString   (a_name);
    // Give name ???
    char*            command;
    command          = CStringCreateF(37+strlen(a_gs),
				      "/vis~/create_view/new_graphics_system %s",
				      a_gs);
    UI->ApplyCommand (command);
    CStringDelete    (command);
    widget           = (Widget)interactorManager->GetCreatedInteractor ();
    if(widget==a_parent) {
      widget = NULL; //No viewer created !
    }
    interactorManager->SetParentInteractor  (NULL);
    interactorManager->SetCreatedInteractor (NULL);
    interactorManager->SetCreationString    (NULL);

  }

  if(widget==NULL) {
    Arg    args[1];
    CWarnF ("G4UIWo::CreateG4View : can't create a %s view.\n Create a label instead.\n",a_gs);
#ifdef HAS_XM
    widget        = XmCreateLabel (a_parent,a_name,args,0);
    XtManageChild (widget);
#endif
#ifdef HAS_XAW
    widget        = XtCreateManagedWidget (a_name,labelWidgetClass,a_parent,args,0);
#endif
  }

  return widget;
}

#include "G4oCommon.icc"

#endif
