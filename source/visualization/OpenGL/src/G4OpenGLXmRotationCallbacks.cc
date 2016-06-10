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
// $Id: G4OpenGLXmRotationCallbacks.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// Andrew Walkden  16th April 1997
// G4OpenGLXmRotationCallbacks : 
//                       Several callback functions used by
//                       elements of the control panel to
//                       rotate the view.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLXmViewer.hh"

#include "G4PhysicalConstants.hh"
#include "G4Scene.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

void G4OpenGLXmViewer::theta_rotation_callback (Widget w, 
					      XtPointer clientData, 
					      XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->rotate_right = get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    rotate_in_theta (pView, NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->rotation_timer);
  }
}

void G4OpenGLXmViewer::rotate_in_theta (XtPointer clientData,
				      XtIntervalId* timer_id) 
{
  //theta spin stuff here
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  if (pView->rotate_right) {
    pView->rotateScene(1,0);
  } else {
    pView->rotateScene(-1,0);
  }
  /*
  G4double delta_theta;

  if (pView->fVP.GetLightsMoveWithCamera()) {
    if (pView->rotate_right) {
      delta_theta = -(pView->fRot_sens);
    } else {
      delta_theta = pView->fRot_sens;
    }
  } else {
    if (pView->rotate_right) {
      delta_theta = pView->fRot_sens;
    } else {
      delta_theta = -(pView->fRot_sens);
    }
  }    
  delta_theta *= deg;
  // Rotates by fixed azimuthal angle delta_theta.

  const G4Vector3D& vp = pView->fVP.GetViewpointDirection ().unit ();
  const G4Vector3D& up = pView->fVP.GetUpVector ().unit ();
  const G4Vector3D& zprime = up;
  G4double cosalpha = up.dot (vp);
  G4double sinalpha = std::sqrt (1. - std::pow (cosalpha, 2));
  G4Vector3D yprime = (zprime.cross (vp)).unit ();
  G4Vector3D xprime = yprime.cross (zprime);
  // Projection of vp on plane perpendicular to up...
  G4Vector3D a1 = sinalpha * xprime;
  // Required new projection...
  G4Vector3D a2 =
    sinalpha * (std::cos (delta_theta) * xprime + std::sin (delta_theta) * yprime);
  // Required Increment vector...
  G4Vector3D delta = a2 - a1;
  // So new viewpoint is...
  G4Vector3D viewPoint = vp + delta;

  pView->fVP.SetViewAndLights (viewPoint);  
  */

  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();

  pView->rotation_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     rotate_in_theta,
     pView);
}  

void G4OpenGLXmViewer::phi_rotation_callback (Widget w, 
					    XtPointer clientData, 
					    XtPointer callData) 
{
  XmArrowButtonCallbackStruct *cbs = (XmArrowButtonCallbackStruct*) callData;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  pView->rotate_up = get_boolean_userData (w);

  if (cbs->reason == XmCR_ARM) {
    rotate_in_phi (pView, NULL);
  } else if (cbs->reason == XmCR_DISARM) {
    XtRemoveTimeOut (pView->rotation_timer);
  }
}

void G4OpenGLXmViewer::rotate_in_phi (XtPointer clientData,
				    XtIntervalId* timer_id)
{
  //phi spin stuff here
  //  G4double delta_alpha;
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*) clientData;

  if (pView -> rotate_up) {
    pView->rotateScene(0,-1);
  } else {
    pView->rotateScene(0,1);
  }
  /*
  if (pView->fVP.GetLightsMoveWithCamera()) {
    if (pView -> rotate_up) {
      delta_alpha = -(pView->fRot_sens);
    } else {
      delta_alpha = pView->fRot_sens;
    }
  } else {
    if (pView -> rotate_up) {
      delta_alpha = pView->fRot_sens;
    } else {
      delta_alpha = -(pView->fRot_sens);
    }
  }    

  delta_alpha *= deg;

  const G4Vector3D& vp = pView->fVP.GetViewpointDirection ().unit ();
  const G4Vector3D& up = pView->fVP.GetUpVector ().unit ();

  const G4Vector3D& xprime = vp;
  G4Vector3D yprime = (up.cross(xprime)).unit();
  G4Vector3D zprime = (xprime.cross(yprime)).unit();

  G4Vector3D new_vp = std::cos(delta_alpha) * xprime + std::sin(delta_alpha) * zprime;

  pView->fVP.SetViewAndLights (new_vp.unit());

  if (pView->fVP.GetLightsMoveWithCamera()) {
    G4Vector3D new_up = (new_vp.cross(yprime)).unit();
    pView->fVP.SetUpVector(new_up);
  }

  */
  pView->SetView ();
  pView->ClearView ();
  pView->DrawView ();
  
  pView->rotation_timer = XtAppAddTimeOut 
    (pView->app,
     timer_id == NULL ? 500 : 1,
     rotate_in_phi,
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
    G4Exception
      ("G4OpenGLXmViewer::set_rot_sens_callback",
       "opengl2004", FatalException,
       "Bad value returned for dp in set_rot_sens_callback");
  }

  pView->fRot_sens = (G4float)(cbs->value) / ten_to_the_dp;
}  

void G4OpenGLXmViewer::set_rot_subject_callback (Widget w, 
					       XtPointer clientData, 
					       XtPointer) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  
  G4int choice = get_int_userData (w);
  
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
      G4Exception
	("G4OpenGLXmViewer::set_rot_subject_callback",
	 "opengl2005", FatalException,
	 "Unrecognised choice made in set_rot_subject_callback"); 
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
     wobble_timer_callback,
     pView);
}  

void G4OpenGLXmViewer::wobble_timer_callback (XtPointer clientData,
					    XtIntervalId*) 
{
  G4OpenGLXmViewer* pView = (G4OpenGLXmViewer*)clientData;
  const G4Vector3D& up = pView->fVP.GetUpVector();
  G4Vector3D third_axis = up.cross(pView->original_vp);
  G4double pi_div_by_ten = pi / 10.0;
  G4Vector3D d_up = 0.1 * (std::sin ((G4double)pView->frameNo * pi_div_by_ten * 2.)) * up;
  G4Vector3D d_third_axis = 0.1 * (std::sin ((G4double)pView->frameNo * (pi_div_by_ten))) * third_axis;

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
       wobble_timer_callback,
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
