//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4WoViewer.cc,v 1.6.4.1 2001/06/28 19:15:36 gunter Exp $
// GEANT4 tag $Name:  $
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
#include "G4Scene.hh"
#include "G4Wo.hh"
#include "G4GoSceneHandler.hh"
//This
#include "G4WoViewer.hh"

/***************************************************************************/
G4WoViewer::G4WoViewer (
 G4GoSceneHandler& scene,
 const G4String& name
)
:G4VViewer   (scene, scene.IncrementViewCount(), name)
,fSceneHandler    (scene)
,fGoCamera (NULL)
,fXoCamera (NULL)
,fShell    (NULL)
/***************************************************************************/
// Constructor. 
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
#ifdef DEBUG
  G4cout << "G4WoViewer::G4WoViewer" << G4endl;
#endif

  if(WoIsInitialized()==0) {
    CWarn ("G4WoViewer : Wo is not initialized !\n");
    return;
  } 
  char* defNameGS = CStringDuplicate((char*)fName.data());
  char* defName = CStringGetFirstWord(defNameGS);
  CStringDelete(defNameGS);
  
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
G4WoViewer::~G4WoViewer (
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
void G4WoViewer::ClearView (
) 
/***************************************************************************/
{
}
/***************************************************************************/
void G4WoViewer::SetView (
) 
/***************************************************************************/
// Calculates view representation based on extent of object being
// viewed and (initial) direction of camera.  (Note: it can change
// later due to user interaction via visualization system's GUI.)
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  const G4Point3D  target
    = fSceneHandler.GetScene()->GetStandardTargetPoint()
    + fVP.GetCurrentTargetPoint ();                                            
  G4double radius  = fSceneHandler.GetScene()->GetExtent().GetExtentRadius();
  if(radius<=0.) radius = 1. * m;
  const G4double cameraDistance = fVP.GetCameraDistance (radius);
  const G4Point3D& pCamera = 
    target + cameraDistance * fVP.GetViewpointDirection().unit();

#ifdef DEBUG
  G4cout << "G4WoViewer::SetView : target.x " << target.x() << " target.y " << target.y() << " target.z " << target.z() << G4endl;
  G4cout << "G4WoViewer::SetView :  dir.x " <<  fVP.GetViewpointDirection().unit().x() << 
                             "  dir.y " <<  fVP.GetViewpointDirection().unit().y() << 
                             "  dir.z " <<  fVP.GetViewpointDirection().unit().z() << G4endl;
  G4cout << "G4WoViewer::SetView : cameraDistance " << cameraDistance << G4endl;
  G4cout << "G4WoViewer::SetView : pCamera.x " << pCamera.x() << " pCamera.y " << pCamera.y() << " pCamera.z " << pCamera.z() << G4endl;
  G4cout << "G4WoViewer::SetView : pCamera.x " << pCamera.x() << " pCamera.y " << pCamera.y() << " pCamera.z " << pCamera.z() << G4endl;
#endif

  OCameraSetCenter          (fGoCamera,target.x() ,target.y() ,target.z());
  OCameraSetDefaultVRP      (fGoCamera,pCamera.x(),pCamera.y(),pCamera.z());

  OCameraSetDefaultUpVector (fGoCamera,0.,1.,0.);

  OCameraSetField           (fGoCamera,-radius,radius);

}
/***************************************************************************/
void G4WoViewer::DrawView (
) 
/***************************************************************************/
{
  NeedKernelVisit ();
  ProcessView     ();
  ShowView        ();
}
/***************************************************************************/
void G4WoViewer::ShowView (
) 
/***************************************************************************/
{
  OCameraViewNode (fGoCamera,fSceneHandler.GetRootNode());  
}
/***************************************************************************/
OCamera G4WoViewer::GetCamera (
) 
/***************************************************************************/
{
  return fGoCamera;
}
/***************************************************************************/
void G4WoViewer::FinishView (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}


#endif
