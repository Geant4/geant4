// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmPanningCallbacks.cc,v 1.1 1999-01-07 16:15:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "G4OpenGLXmView.hh"

void G4OpenGLXmView::zoom_callback (Widget w, 
				    XtPointer clientData, 
				    XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;
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
    G4cout << "dp is " << dp << endl;
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

void G4OpenGLXmView::dolly_callback (Widget w, 
				     XtPointer clientData, 
				     XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;
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
    G4cout << "dp is " << dp << endl;
    return;
  }

  G4double dolly = (G4double)(cbs->value) / ten_to_the_dp;

  pView->fVP.IncrementDolly (dolly);
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

}  

void G4OpenGLXmView::pan_left_right_callback (Widget w, 
					      XtPointer clientData, 
					      XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;

  pView->pan_right = G4OpenGLXmView::get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    left_right_pan_callback (pView,NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->pan_timer);
  }
}

void G4OpenGLXmView::left_right_pan_callback (XtPointer clientData,
					      XtIntervalId* timer_id) 

{
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;
  G4double delta;
  
  if (pView->pan_right) {
    delta = (G4double)pView->pan_sens;
  } else {
    delta = -((G4double)pView->pan_sens);
  }
  
  G4Point3D tp = pView -> fVP.GetCurrentTargetPoint ();
  
  const G4Vector3D& upVector = pView->fVP.GetUpVector ();
  const G4Vector3D& vpVector = pView->fVP.GetViewpointDirection ();
  
  G4Vector3D unitRight = (upVector.cross (vpVector)).unit();
  G4Vector3D unitUp    = (vpVector.cross (unitRight)).unit();
  
  tp += delta * unitRight;
  pView->fVP.SetCurrentTargetPoint (tp);

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->pan_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     left_right_pan_callback,
     pView);
}  

void G4OpenGLXmView::pan_up_down_callback (Widget w, 
					   XtPointer clientData, 
					   XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;

  pView->pan_up = G4OpenGLXmView::get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    up_down_pan_callback (pView,NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->pan_timer);
  }
}

void G4OpenGLXmView::up_down_pan_callback (XtPointer clientData,
					   XtIntervalId* timer_id) 
{
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;
  G4double delta;
  
  if (pView->pan_up) {
    delta = (G4double)pView->pan_sens;
  } else {
    delta = -((G4double)pView->pan_sens);
  }
  
  G4Point3D tp = pView->fVP.GetCurrentTargetPoint ();
  const G4Vector3D& upVector = pView->fVP.GetUpVector ();
  const G4Vector3D& vpVector = pView->fVP.GetViewpointDirection ();
  
  G4Vector3D unitRight = (upVector.cross (vpVector)).unit();
  G4Vector3D unitUp    = (vpVector.cross (unitRight)).unit();
  tp += delta * unitUp;
  pView->fVP.SetCurrentTargetPoint (tp);
  
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->pan_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     up_down_pan_callback,
     pView);
}  

void G4OpenGLXmView::set_pan_sens_callback (Widget w, 
					    XtPointer clientData, 
					    XtPointer callData) 
{
  XmScaleCallbackStruct *cbs = (XmScaleCallbackStruct*) callData;
  G4OpenGLXmView* pView = (G4OpenGLXmView*) clientData;
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
    G4cout << "dp is " << dp << endl;
    return;
  }
  
  pView->pan_sens = (G4double)((cbs->value) / ten_to_the_dp);
}  

#endif


