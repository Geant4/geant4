// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmRotationCallbacks.hh,v 1.1 1999-01-07 16:14:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMROTATIONCALLBACKS_HH
#define G4OPENGLXMROTATIONCALLBACKS_HH

#include "G4OpenGLXmView.hh"

void theta_rotation_callback (Widget w, 
			      XtPointer clientData, 
			      XtPointer callData);

void rotate_in_theta (XtPointer clientData, 
		      XtIntervalId timer_id);

void phi_rotation_callback (Widget w, 
			    XtPointer clientData, 
			    XtPointer callData);

void rotate_in_phi (XtPointer clientData, 
		    XtIntervalId timer_id);

void set_rot_sens_callback (Widget w, 
			    XtPointer clientData, 
			    XtPointer callData);

void set_rot_subject_callback (Widget w, 
			       XtPointer clientData, 
			       XtPointer callData);

#endif

#endif
