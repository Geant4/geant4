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
// $Id: G4OpenGLXmPanningCallbacks.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// Andrew Walkden  16th April 1997
// G4OpenGLXmPanningCallbacks : 
//                       Several callback functions used by
//                       elements of the control panel to`pan' 
//                       the view (i.e. move the viewpoint and 
//                       camera positions by equal amounts).
//                       Zoom callback is also here.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"
#include "G4VSceneHandler.hh"
#include <Xm/ToggleB.h>

#include "G4Scene.hh"

void G4OpenGLXmViewer::zoom_callback (Widget w, 
				    XtPointer clientData, 
				    XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  short dp = -1;
  G4float ten_to_the_dp = 10.;

  XtVaGetValues (w, 
		 XmNdecimalPoints, &dp,
		 NULL);

  if (dp == 0) {
    ten_to_the_dp = 1.;
  } else if ( dp > 0) {
    for (G4int i = 1; i < (G4int)dp; i++) {
      ten_to_the_dp *= 10.;
    }
  } else {
    G4cout << "dp is " << dp << G4endl;
    return;
  }


  G4double zoomBy = (G4double)(cbs->value) / ten_to_the_dp;
  if (zoomBy <= 0.01) {
    zoomBy = 0.01;
  }

  pView->fVP.SetZoomFactor (zoomBy);
  pView->SetView ();
  pView->ClearView ();
  pView -> DrawView ();
}  

void G4OpenGLXmViewer::dolly_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  short dp = -1;
  G4float ten_to_the_dp = 10.;

  XtVaGetValues (w, 
		 XmNdecimalPoints, &dp,
		 NULL);

  if (dp == 0) {
    ten_to_the_dp = 1.;
  } else if ( dp > 0) {
    for (G4int i = 1; i < (G4int)dp; i++) {
      ten_to_the_dp *= 10.;
    }
  } else {
    G4cout << "dp is " << dp << G4endl;
    return;
  }

  G4double dolly = (G4double)(cbs->value) / ten_to_the_dp;

  pView->fVP.SetDolly (dolly);
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

}  

void G4OpenGLXmViewer::pan_left_right_callback (Widget w, 
					      XtPointer clientData, 
					      XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->pan_right = get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    left_right_pan_callback (pView,NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->pan_timer);
  }
}

void G4OpenGLXmViewer::left_right_pan_callback (XtPointer clientData,
					      XtIntervalId* timer_id) 

{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  G4double delta;
  
  if (pView->pan_right) {
    delta = pView->fPan_sens;
  } else {
    delta = -pView->fPan_sens;
  }
  
  G4Point3D stp
    = pView -> GetSceneHandler()->GetScene()->GetStandardTargetPoint();

  G4Point3D tp = stp + pView -> fVP.GetCurrentTargetPoint ();
  
  const G4Vector3D& upVector = pView->fVP.GetUpVector ();
  const G4Vector3D& vpVector = pView->fVP.GetViewpointDirection ();
  
  G4Vector3D unitRight = (upVector.cross (vpVector)).unit();
  G4Vector3D unitUp    = (vpVector.cross (unitRight)).unit();
  
  tp += delta * unitRight;
  pView->fVP.SetCurrentTargetPoint (tp - stp);

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->pan_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     left_right_pan_callback,
     pView);
}  

void G4OpenGLXmViewer::pan_up_down_callback (Widget w, 
					   XtPointer clientData, 
					   XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->pan_up = get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    up_down_pan_callback (pView,NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->pan_timer);
  }
}

void G4OpenGLXmViewer::up_down_pan_callback (XtPointer clientData,
					   XtIntervalId* timer_id) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  G4double delta;
  
  if (pView->pan_up) {
    delta = pView->fPan_sens;
  } else {
    delta = -(pView->fPan_sens);
  }
  
  G4Point3D stp
    = pView -> GetSceneHandler()->GetScene()->GetStandardTargetPoint();
  G4Point3D tp = stp + pView -> fVP.GetCurrentTargetPoint ();
  const G4Vector3D& upVector = pView->fVP.GetUpVector ();
  const G4Vector3D& vpVector = pView->fVP.GetViewpointDirection ();
  
  G4Vector3D unitRight = (upVector.cross (vpVector)).unit();
  G4Vector3D unitUp    = (vpVector.cross (unitRight)).unit();
  tp += delta * unitUp;
  pView->fVP.SetCurrentTargetPoint (tp - stp);
  
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->pan_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     up_down_pan_callback,
     pView);
}  

void G4OpenGLXmViewer::set_pan_sens_callback (Widget w, 
					    XtPointer clientData, 
					    XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  short dp = -1;
  G4float ten_to_the_dp = 10.;

  XtVaGetValues (w, 
		 XmNdecimalPoints, &dp,
		 NULL);

  if (dp == 0) {
    ten_to_the_dp = 1.;
  } else if ( dp > 0) {
    for (G4int i = 1; i < (G4int)dp; i++) {
      ten_to_the_dp *= 10.;
    }
  } else {
    G4cout << "dp is " << dp << G4endl;
    return;
  }
  
  pView->fPan_sens = (G4double)((cbs->value) / ten_to_the_dp);
}  

#endif


