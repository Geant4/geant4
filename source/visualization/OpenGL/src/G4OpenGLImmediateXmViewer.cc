// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXmViewer.cc,v 1.1 1999-01-09 16:23:15 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// Class G4OpenGLImmediateXmViewer : a class derived from G4OpenGLXmViewer
//                                     and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLImmediateXmViewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

G4OpenGLImmediateXmViewer::G4OpenGLImmediateXmViewer
(G4OpenGLImmediateSceneHandler& scene,
 const G4String& name):
G4OpenGLViewer (scene),
G4OpenGLXmViewer (scene),
G4OpenGLImmediateViewer (scene),
G4VViewer (scene, scene.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.

// ensure a suitable window was found
  if (!vi_immediate) {
    G4cerr << "G4OpenGLImmediateXmViewer::G4OpenGLImmediateXmViewer -"
      " G4OpenGLXmViewer couldn't get a visual." << endl;  
    fViewId = -1;  // This flags an error.
    return;
  }

  CreateGLXContext (vi_immediate);

  InitializeGLView ();

  CreateMainWindow ();

  // clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
  
  glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glLineWidth (1.0);

  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glShadeModel (GL_FLAT);

  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //The following code is useless in its current position, as the 
  //G4OpenGLXmViewer constructor gets called *after* it, and hence sets
  //doublebuffer to true or false there, after our little test to correct
  //it in the case of a double buffer being got for an immediate view.
  //Hence, code moved to DrawView.
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  //  if (doublebuffer) {
  //    doublebuffer = false;
  //    glDrawBuffer (GL_FRONT);
  //  }
}

void G4OpenGLImmediateXmViewer::DrawView () {

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
