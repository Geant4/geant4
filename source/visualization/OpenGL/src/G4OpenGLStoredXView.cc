// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredXView.cc,v 1.1 1999-01-07 16:14:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredXView : a class derived from G4OpenGLXView and
//                             G4OpenGLStoredView.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLStoredXView.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

#include <X11/Xatom.h>
#include <X11/Xutil.h>

G4OpenGLStoredXView::G4OpenGLStoredXView (G4OpenGLStoredScene& scene,
					  const G4String& name):
G4OpenGLView (scene),
G4OpenGLXView (scene),
G4OpenGLStoredView (scene),
G4VView (scene, scene.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.

  if (!vi_stored) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLStoredXView::G4OpenGLStoredXView -"
      " G4OpenGLXView couldn't get a visual." << endl;
    return;
  }

  CreateGLXContext (vi_stored);

  InitializeGLView ();

  CreateMainWindow ();

// clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);

}

void G4OpenGLStoredXView::DrawView () {

  if (white_background == true) {
    glClearColor (1., 1., 1., 1.);
  } else {
    glClearColor (0., 0., 0., 1.);
  }

  //Make sure current viewer is attached and clean...
  glXMakeCurrent (dpy, win, cx);
  glViewport (0, 0, WinSize_x, WinSize_y);
  ClearView ();

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  //See if things have changed from last time and remake if necessary...
  KernelVisitDecision ();
  ProcessView ();

  if(style!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    DrawDisplayLists ();
    glFlush ();

    HaloingSecondPass ();

  }

  DrawDisplayLists ();

  FinishView ();
  
}

#endif
