// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmRotationCallbacks.cc,v 1.7 2000-05-13 10:48:03 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  16th April 1997
// G4OpenGLXmRotationCallbacks : 
//                       Several callback functions used by
//                       elements of the control panel to
//                       rotate the view.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"

#include "G4Scene.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

void G4OpenGLXmViewer::theta_rotation_callback (Widget w, 
					      XtPointer clientData, 
					      XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->rotate_right = G4OpenGLXmViewer::get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    G4OpenGLXmViewer::rotate_in_theta (pView, NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->rotation_timer);
  }
}

void G4OpenGLXmViewer::rotate_in_theta (XtPointer clientData,
				      XtIntervalId* timer_id) 
{
  //theta spin stuff here
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;
  G4double delta_theta;

  if (pView->fVP.GetLightsMoveWithCamera()) {
    if (pView->rotate_right) {
      delta_theta = -((G4double)pView->rot_sens);
    } else {
      delta_theta = (G4double)pView->rot_sens;
    }
  } else {
    if (pView->rotate_right) {
      delta_theta = (G4double)pView->rot_sens;
    } else {
      delta_theta = -((G4double)pView->rot_sens);
    }
  }    
  delta_theta *= deg;
  // Rotates by fixed azimuthal angle delta_theta.

  G4Vector3D vp = pView->fVP.GetViewpointDirection ().unit ();
  G4Vector3D up = pView->fVP.GetUpVector ().unit ();
  G4Vector3D& zprime = up;
  G4double cosalpha = up.dot (vp);
  G4double sinalpha = sqrt (1. - pow (cosalpha, 2));
  G4Vector3D viewPoint = pView->fVP.GetViewpointDirection ().unit ();
  G4Vector3D yprime = (zprime.cross (viewPoint)).unit ();
  G4Vector3D xprime = yprime.cross (zprime);
  // Projection of vp on plane perpendicular to up...
  G4Vector3D a1 = sinalpha * xprime;
  // Required new projection...
  G4Vector3D a2 =
    sinalpha * (cos (delta_theta) * xprime + sin (delta_theta) * yprime);
  // Required Increment vector...
  G4Vector3D delta = a2 - a1;
  // So new viewpoint is...
  viewPoint += delta;

  pView->fVP.SetViewAndLights (viewPoint);  

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->rotation_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     G4OpenGLXmViewer::rotate_in_theta,
     pView);
}  

void G4OpenGLXmViewer::phi_rotation_callback (Widget w, 
					    XtPointer clientData, 
					    XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->rotate_up = G4OpenGLXmViewer::get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    G4OpenGLXmViewer::rotate_in_phi (pView, NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->rotation_timer);
  }
}

void G4OpenGLXmViewer::rotate_in_phi (XtPointer clientData,
				    XtIntervalId* timer_id)
{
  //phi spin stuff here
  G4double delta_alpha;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  if (pView->fVP.GetLightsMoveWithCamera()) {
    if (pView -> rotate_up) {
      delta_alpha = -((G4double)pView->rot_sens);
    } else {
      delta_alpha = (G4double)pView->rot_sens;
    }
  } else {
    if (pView -> rotate_up) {
      delta_alpha = (G4double)pView->rot_sens;
    } else {
      delta_alpha = -((G4double)pView->rot_sens);
    }
  }    

  delta_alpha *= deg;

  G4Vector3D vp = pView->fVP.GetViewpointDirection ().unit ();
  G4Vector3D up = pView->fVP.GetUpVector ().unit ();

  G4Vector3D xprime = vp;
  G4Vector3D yprime = (up.cross(xprime)).unit();
  G4Vector3D zprime = (xprime.cross(yprime)).unit();

  G4Vector3D new_vp = cos(delta_alpha) * xprime + sin(delta_alpha) * zprime;

  pView->fVP.SetViewAndLights (new_vp.unit());

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
  
  pView->rotation_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     G4OpenGLXmViewer::rotate_in_phi,
     pView);
}

void G4OpenGLXmViewer::set_rot_sens_callback (Widget w, 
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
    G4Exception("Bad value returned for dp in set_rot_sens_callback");
  }

  pView->rot_sens = (G4float)(cbs->value) / ten_to_the_dp;
}  

void G4OpenGLXmViewer::set_rot_subject_callback (Widget w, 
					       XtPointer clientData, 
					       XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = G4OpenGLXmViewer::get_int_userData (w);
  
  switch (choice) {
  case 0: 
    {
      pView->fVP.SetLightsMoveWithCamera (true);
      break;
    }
  case 1:
    {
      pView->fVP.SetLightsMoveWithCamera (false);
      break;
    }
  default:
    {
      G4Exception("Unrecognised choice made in set_rot_subject_callback"); 
    }
  }
}  

void G4OpenGLXmViewer::wobble_callback (Widget w, 
				      XtPointer, 
				      XtPointer) 
{
  G4OpenGLXmViewer* pView;
  
  XtVaGetValues (w,
		 XmNuserData, &pView,
		 NULL);
  
  pView->original_vp = pView->fVP.GetViewpointDirection();
  pView->wobble_timer = XtAppAddTimeOut
    (pView->app,
     (long unsigned int) (1000. * (1. / pView->wob_sens)),
     G4OpenGLXmViewer::wobble_timer_callback,
     pView);
}  

void G4OpenGLXmViewer::wobble_timer_callback (XtPointer clientData,
					    XtIntervalId*) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  G4Vector3D up = pView->fVP.GetUpVector();
  G4Vector3D vp = pView->fVP.GetViewpointDirection();
  G4Vector3D third_axis = up.cross(pView->original_vp);
  G4double pi_div_by_ten = M_PI / 10.0;
  G4Vector3D d_up = 0.1 * (sin ((G4double)pView->frameNo * pi_div_by_ten * 2.)) * up;
  G4Vector3D d_third_axis = 0.1 * (sin ((G4double)pView->frameNo * (pi_div_by_ten))) * third_axis;

  pView->fVP.SetViewAndLights (pView->original_vp + d_up + d_third_axis);

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
  
  if (pView->frameNo++ == 20) {
    if (pView->wobble_timer) {
      XtRemoveTimeOut (pView->wobble_timer);
      pView->frameNo = 0;
      pView->fVP.SetViewAndLights (pView->original_vp);
      pView->SetView ();
      pView->ClearView ();
      pView->DrawView ();
   }
  } else {
    pView->wobble_timer = XtAppAddTimeOut
      (pView->app,
       (long unsigned int) (1000. * (1. / pView->wob_sens)),
       G4OpenGLXmViewer::wobble_timer_callback,
       pView);
  }
}

void G4OpenGLXmViewer::reset_callback (Widget w, 
				     XtPointer, 
				     XtPointer) 
{
  
  G4OpenGLXmViewer* pView;
  
  XtVaGetValues (w,
		 XmNuserData, &pView,
		 NULL);
  
  pView->fVP.SetCurrentTargetPoint(G4Point3D());
  pView->fVP.SetZoomFactor(1.0);
  pView->fVP.SetDolly(0.0);
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
  pView->zoom_low = 0.1;
  pView->zoom_high = 10.0;
  
}
#endif
