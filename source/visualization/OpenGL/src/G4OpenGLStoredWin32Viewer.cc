// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredWin32Viewer.cc,v 1.1 1999-01-09 16:23:21 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLStoredWin32Viewer : a class derived from G4OpenGLWin32Viewer and
//                             G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLStoredWin32Viewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

G4OpenGLStoredWin32Viewer::G4OpenGLStoredWin32Viewer (G4OpenGLStoredSceneHandler& scene):
G4OpenGLViewer (scene),
G4OpenGLWin32Viewer (scene),
G4OpenGLStoredViewer (scene),
G4VViewer (scene, scene.IncrementViewCount ()) {

  if (fViewId < 0) return;  // In case error in base class instantiation.

  //Check that G4OpenGLWin32Viewer got a double buffered colour visual

  CreateGLWin32Context ();
  CreateMainWindow ();

// clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);

}

void G4OpenGLStoredWin32Viewer::DrawView () {

  if (white_background == true) {
    glClearColor (1., 1., 1., 1.);
  } else {
    glClearColor (0., 0., 0., 1.);
  }

  //Make sure current viewer is attached and clean...
  //Win32 version needed
  //  glXMakeCurrent (dpy, win, cx);
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
