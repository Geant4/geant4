// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateWin32View.cc,v 1.1 1999-01-07 16:14:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLImmediateWin32View : a class derived from G4OpenGLWin32View and
//                                G4OpenGLImmediateView.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLImmediateWin32View.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

G4OpenGLImmediateWin32View::G4OpenGLImmediateWin32View (G4OpenGLImmediateScene& scene):
G4OpenGLView (scene),
G4OpenGLWin32View (scene),
G4OpenGLImmediateView (scene),
G4VView (scene, scene.IncrementViewCount ()) {

  if (fViewId < 0) return;  // In case error in base class instantiation.

// ensure a suitable window was found

  CreateGLWin32Context ();
  CreateMainWindow ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  if (doublebuffer) {
    doublebuffer = false;
    glDrawBuffer (GL_FRONT);
  }

  // clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
  
}

void G4OpenGLImmediateWin32View::DrawView () {

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
  //Win32 version needed
  //glXMakeCurrent (dpy, win, cx);
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
