// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4WoView.cc,v 2.8 1998/11/09 15:51:38 barrand Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Guy Barrand 04 November 1996

#ifdef G4VIS_BUILD_OPACS_DRIVER

//#define DEBUG

/*Co*/
#include <CPrinter.h>
#include <CString.h>
/*Xx*/
#include <XWidget.h>
/*Go*/
#include <OCamera.h>
/*Wo*/
#include <OWidget.h>
#include <Wo.h>
//G4
#include "G4VInteractorManager.hh"
#include "G4Wo.hh"
#include "G4GoScene.hh"
//This
#include "G4WoView.hh"

/***************************************************************************/
G4WoView::G4WoView (
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
  G4cout << "G4WoView::G4WoView" << endl;
#endif

  if(WoIsInitialized()==0) {
    CWarn ("G4WoView : Wo is not initialized !\n");
    return;
  } 

  char* defName = CStringDuplicate((char*)fName.data());
  CStringReplacePart (&defName,"-","_");

  fGoCamera = OCameraGetIdentifier (defName);
  if(fGoCamera!=NULL) {
    return;  //OCamera could have been created 
             //by loading an odb file at Wo startup.
  }

  G4VInteractorManager* interactorManager = G4Wo::GetInteractorManager();
  Widget parent = (Widget)interactorManager->GetParentInteractor ();
  char* wname = NULL;
  if(parent==NULL) { //Create a shell widget.
    wname = CStringCreateF (6+strlen(defName),"shell_%s",defName);
#ifdef HAS_XM
    fShell = OWidgetCreate  (wname,XWidgetGetTop(),"XmFormDialog",False);
#else
    fShell = OWidgetCreate  (wname,XWidgetGetTop(),"TopLevelShell",False);
#endif
    CStringDelete (wname);   
    parent        = fShell;
    wname          = CStringDuplicate (defName);
  } else {
    char* str = interactorManager->GetCreationString ();
    if(str==NULL) {
      wname = CStringDuplicate (defName);
    } else {
      wname = CStringDuplicate (str);
    }
  }


  fXoCamera = OWidgetCreate  (wname,parent,"XoCamera",False);

  if(fShell!=NULL) {
#ifdef HAS_XM
    OWidgetSetResourceFromString (fShell,"dialogTitle",wname,False);
    OWidgetSetResourceFromString (fXoCamera,"topAttachment","attach_form",False);
    OWidgetSetResourceFromString (fXoCamera,"leftAttachment","attach_form",False);
    OWidgetSetResourceFromString (fXoCamera,"rightAttachment","attach_form",False);
    OWidgetSetResourceFromString (fXoCamera,"bottomAttachment","attach_form",False);
#else
    OWidgetSetResourceFromString (fShell,"title",wname,False);
#endif
    XWidgetMap    (fShell);
  }

  fGoCamera     = OCameraGetIdentifier (wname);

  if(fGoCamera==NULL) {
    CWarnF ("Can't find/create OCamera %s\n",wname);
  }

  CStringDelete (wname);      
  CStringDelete (defName);      

  interactorManager->SetCreatedInteractor (fXoCamera);
}
/***************************************************************************/
G4WoView::~G4WoView (
) 
/***************************************************************************/
// Destructor.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(fShell!=NULL) {
    OWidgetDelete (fShell,False);
  } else if(fXoCamera!=NULL) {
    OWidgetDelete (fXoCamera,False);
  }
}
/***************************************************************************/
void G4WoView::ClearView (
) 
/***************************************************************************/
{
}
/***************************************************************************/
void G4WoView::SetView (
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
  G4cout << "G4WoView::SetView : target.x " << target.x() << " target.y " << target.y() << " target.z " << target.z() << endl;
  G4cout << "G4WoView::SetView :  dir.x " <<  fVP.GetViewpointDirection().unit().x() << 
                             "  dir.y " <<  fVP.GetViewpointDirection().unit().y() << 
                             "  dir.z " <<  fVP.GetViewpointDirection().unit().z() << endl;
  G4cout << "G4WoView::SetView : cameraDistance " << cameraDistance << endl;
  G4cout << "G4WoView::SetView : pCamera.x " << pCamera.x() << " pCamera.y " << pCamera.y() << " pCamera.z " << pCamera.z() << endl;
  G4cout << "G4WoView::SetView : pCamera.x " << pCamera.x() << " pCamera.y " << pCamera.y() << " pCamera.z " << pCamera.z() << endl;
#endif

  OCameraSetCenter          (fGoCamera,target.x() ,target.y() ,target.z());
  OCameraSetDefaultVRP      (fGoCamera,pCamera.x(),pCamera.y(),pCamera.z());

  OCameraSetDefaultUpVector (fGoCamera,0.,1.,0.);

  OCameraSetField           (fGoCamera,-radius,radius);

}
/***************************************************************************/
void G4WoView::DrawView (
) 
/***************************************************************************/
{
  NeedKernelVisit ();
  ProcessView     ();
  ShowView        ();
}
/***************************************************************************/
void G4WoView::ShowView (
) 
/***************************************************************************/
{
  OCameraViewNode (fGoCamera,fScene.GetRootNode());  
}
/***************************************************************************/
OCamera G4WoView::GetCamera (
) 
/***************************************************************************/
{
  return fGoCamera;
}
/***************************************************************************/
void G4WoView::FinishView (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}


#endif
