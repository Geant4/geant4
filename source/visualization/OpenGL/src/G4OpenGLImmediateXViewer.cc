// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXViewer.cc,v 1.3 1999-12-15 14:54:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateXViewer : a class derived from G4OpenGLXViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLImmediateXViewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

#include <X11/Xatom.h>
#include <X11/Xutil.h>

G4OpenGLImmediateXViewer::G4OpenGLImmediateXViewer (G4OpenGLImmediateSceneHandler& scene,
						const G4String&  name):
G4OpenGLViewer (scene),
G4OpenGLXViewer (scene),
G4OpenGLImmediateViewer (scene),
G4VViewer (scene, scene.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.

// ensure a suitable window was found
  if (!vi_immediate) {
    G4cerr << "G4OpenGLImmediateXViewer::G4OpenGLImmediateXViewer -"
      " G4OpenGLXViewer couldn't get a visual." << G4endl;  
    fViewId = -1;  // This flags an error.
    return;
  }

  CreateGLXContext (vi_immediate);

  InitializeGLView ();

  CreateMainWindow ();

  // clear the buffers and window.
  ClearView ();
  FinishView ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  if (doublebuffer) {
    doublebuffer = false;
    glDrawBuffer (GL_FRONT);
  }

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
  
}

G4OpenGLImmediateXViewer::~G4OpenGLImmediateXViewer () {}

void G4OpenGLImmediateXViewer::DrawView () {

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  if (doublebuffer) {
    doublebuffer = false;
    glDrawBuffer (GL_FRONT);
  }

  if (white_background == true) {
    glClearColor (1., 1., 1., 1.);
  } else {
    glClearColor (0., 0., 0., 1.);
  }
  glClearDepth (1.0);

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  //Make sure current viewer is attached and clean...
  glXMakeCurrent (dpy, win, cx);
  glViewport (0, 0, WinSize_x, WinSize_y);
  ClearView ();

  if(style!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    NeedKernelVisit ();
    ProcessView ();
    glFlush ();

    HaloingSecondPass ();

  }

  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
  FinishView ();

}

#endif
