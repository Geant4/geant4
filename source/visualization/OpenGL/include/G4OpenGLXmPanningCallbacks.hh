// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmPanningCallbacks.hh,v 1.1 1999-01-07 16:14:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMPANNINGCALLBACKS_HH
#define G4OPENGLXMPANNINGCALLBACKS_HH

#include "G4OpenGLXmView.hh"

void zoom_callback (Widget w, 
		    XtPointer clientData, 
		    XtPointer callData);

void dolly_callback (Widget w, 
		     XtPointer clientData, 
		     XtPointer callData);

void pan_left_right_callback (Widget w, 
			      XtPointer clientData, 
			      XtPointer callData);

void pan_up_down_callback (Widget w, 
			   XtPointer clientData, 
			   XtPointer callData);

void left_right_pan_callback (XtPointer clientData, 
			      XtIntervalId timer_id);

void up_down_pan_callback (XtPointer clientData, 
			   XtIntervalId timer_id);

#endif

#endif
