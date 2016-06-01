// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4XoView.cc,v 2.7 1998/11/06 13:42:11 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Guy Barrand 04 November 1996

#ifdef G4VIS_BUILD_OPACS_DRIVER

//#define DEBUG

//X11
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
/*Co*/
#include <CPrinter.h>
#include <CString.h>
#include <CList.h>
/*Xx*/
#include <XWidget.h>
/*Go*/
#include <OCamera.h>
/*Xo*/
#include <XoCamera.h>
//G4
#include "G4GoScene.hh"
#include "G4Xo.hh"
#include "G4VInteractorManager.hh"
//This
#include "G4XoView.hh"

static Bool   WaitForNotify    (Display*,XEvent*,char*);
static void   ActivateCallback (Widget,XtPointer,XtPointer);
static void   CollectCallback  (Widget,XtPointer,XtPointer);
/***************************************************************************/
G4XoView::G4XoView (
 G4GoScene& scene,
 const G4String& name
)
:G4VView   (scene, scene.IncrementViewCount(), name)
,fScene    (scene)
,fGoCamera (NULL)
,fXoCamera (NULL)
,fShell    (NULL)
/***************************************************************************/
// Constructor. 
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4XoView::G4XoView" << endl;
#endif

  G4VInteractorManager* interactorManager = G4Xo::GetInteractorManager();
  Widget top    = (Widget)interactorManager->GetMainInteractor ();
  if(top==NULL) {
    CWarn ("Xt is not initialized !\n");
    return;
  } 

  char* defName = CStringDuplicate((char*)fName.data());
  CStringReplacePart (&defName,"-","_");

  Arg    args[1];
  char*  wname = NULL;
  Widget parent = (Widget)interactorManager->GetParentInteractor ();
  if(parent==NULL) {  
    // Create a shell widget.
    // Add WMClose does not works on ApplicationShell, then use a TopLevelShell.
    char*  cname;
    wname  = CStringCreateF   (6+strlen(defName),"shell_%s",defName);
    cname  = CStringCreateF   (6+strlen(defName),"Shell_%s",defName);
    fShell = XtAppCreateShell (wname,cname,topLevelShellWidgetClass,XtDisplay(top),args,0); 
    CStringDelete   (wname);
    CStringDelete   (cname);
    interactorManager->AddShell   (fShell);
    parent          = fShell;
    wname           = CStringDuplicate(defName);
  } else {
    char* str = interactorManager->GetCreationString ();
    if(str==NULL) {
      wname = CStringDuplicate (defName);
    } else {
      wname = CStringDuplicate (str);
    }
  }

  fXoCamera       = XtCreateManagedWidget (wname,xoCameraWidgetClass ,parent, args, 0);
  fGoCamera       = (OCamera) XoCameraGetCamera (fXoCamera);
  
  XoCameraAddPopupEntry (fXoCamera,"Escape",ActivateCallback,NULL);
  XtAddCallback         (fXoCamera,XoNcollectCallback,CollectCallback,NULL);

  if(fShell!=NULL) {
    XtRealizeWidget (fShell);
    XtMapWidget     (fShell);
    XEvent          event;
    XIfEvent        (XtDisplay(fShell), &event, WaitForNotify, (char*)XtWindow(fShell));
  }  
  
  if(fGoCamera==NULL) {
    CWarnF    ("Can't create OCamera %s\n",wname);
    fXoCamera = NULL;
    fGoCamera = NULL;
  }      

  // g++ wants the cast !!!
  interactorManager->AddSecondaryLoopPreAction  ((G4SecondaryLoopAction)XoCameraEnablePopup);
  interactorManager->AddSecondaryLoopPostAction ((G4SecondaryLoopAction)XoCameraDisablePopup);
  interactorManager->SetCreatedInteractor       (fXoCamera);

  CStringDelete   (wname);
  CStringDelete   (defName);
}
/***************************************************************************/
G4XoView::~G4XoView (
) 
/***************************************************************************/
// Destructor.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(fShell!=NULL) {
    G4Xo::GetInteractorManager()->RemoveShell (fShell);
    XtDestroyWidget (fShell);
  } else if(fXoCamera!=NULL) {
    XtDestroyWidget (fXoCamera);
  }
  fXoCamera         = NULL;
  fGoCamera         = NULL;
}
/***************************************************************************/
void G4XoView::ClearView (
) 
/***************************************************************************/
{
}
/***************************************************************************/
void G4XoView::SetView (
) 
/***************************************************************************/
// Calculates view representation based on extent of object being
// viewed and (initial) direction of camera.  (Note: it can change
// later due to user interaction via visualization system's GUI.)
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  const G4Point3D  target  = fVP.GetCurrentTargetPoint ();
  G4double radius  = fScene.GetSceneData().GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1. * m;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Point3D& pCamera = 
    target + cameraDistance * fVP.GetViewpointDirection().unit();

#ifdef DEBUG
  G4cout << "G4XoView::SetView : target.x " << target.x() << " target.y " << target.y() << " target.z " << target.z() << endl;
  G4cout << "G4XoView::SetView :  dir.x " <<  fVP.GetViewpointDirection().unit().x() << 
                             "  dir.y " <<  fVP.GetViewpointDirection().unit().y() << 
                             "  dir.z " <<  fVP.GetViewpointDirection().unit().z() << endl;
  G4cout << "G4XoView::SetView : cameraDistance " << cameraDistance << endl;
  G4cout << "G4XoView::SetView : pCamera.x " << pCamera.x() << " pCamera.y " << pCamera.y() << " pCamera.z " << pCamera.z() << endl;
#endif

  OCameraSetCenter          (fGoCamera,target.x() ,target.y() ,target.z());
  OCameraSetDefaultVRP      (fGoCamera,pCamera.x(),pCamera.y(),pCamera.z());

  OCameraSetDefaultUpVector (fGoCamera,0.,1.,0.);

  OCameraSetField           (fGoCamera,-radius,radius);

}
/***************************************************************************/
OCamera G4XoView::GetCamera (
) 
/***************************************************************************/
{
  return fGoCamera;
}
/***************************************************************************/
void G4XoView::DrawView (
) 
/***************************************************************************/
{
  NeedKernelVisit ();
  ProcessView     ();
  FinishView      ();
}
/***************************************************************************/
void G4XoView::FinishView (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
void G4XoView::ShowView (
) 
/***************************************************************************/
{
  OCameraViewNode (fGoCamera,fScene.GetRootNode());  
  G4Xo::GetInteractorManager() -> SecondaryLoop();
}
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
Bool WaitForNotify (
 Display* d
,XEvent*  e
,char*    arg
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
/*.........................................................................*/
  d = NULL;
  return (e->type == MapNotify) && (e->xmap.window == (Window)arg);
}
/***************************************************************************/
void ActivateCallback (
 Widget    This 
,XtPointer a_tag
,XtPointer a_data
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4Xo::GetInteractorManager() -> RequireExitSecondaryLoop (XO_EXIT_CODE);
  This   = NULL;
  a_tag  = NULL;
  a_data = NULL;
}
/***************************************************************************/
void CollectCallback (
 Widget    This 
,XtPointer a_tag
,XtPointer a_data
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  OCamera camera;
  ONode*  nodes;
  int     noden,count;
/*.........................................................................*/
  // collectCallback = osh> getCameraNodes `thisCamera` | collect ONode where highlight eq 1 | get - name | dump -
  camera = (OCamera)XoCameraGetCamera (This);
  nodes  = OCameraCollect (camera,OCollectHighlighted);
  noden  = CListGetSize   ((CList)nodes);
  for(count=0;count<noden;count++)
    {
      char* name = ONodeGetName(nodes[count]);
      CInfoF("%s\n",name==NULL ? "(nil)" : name);
    }
  CListDelete ((CList)nodes);
  This   = NULL;
  a_tag  = NULL;
  a_data = NULL;
}

#endif
