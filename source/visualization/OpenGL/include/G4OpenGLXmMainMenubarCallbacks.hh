// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLXmMainMenubarCallbacks.hh,v 1.3 1999-12-15 14:54:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMMAINMENUBARCALLBACKS_HH
#define G4OPENGLXMMAINMENUBARCALLBACKS_HH

#include "G4OpenGLXmViewer.hh"

void misc_callback (Widget w, XtPointer clientData, XtPointer callData);
void actions_callback (Widget w, XtPointer clientData, XtPointer callData);
void update_panels_callback (Widget w, XtPointer clientData, XtPointer callData);

#endif

#endif
