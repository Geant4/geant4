//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLImmediateWin32Viewer.cc,v 1.6 2002-10-16 10:44:15 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLImmediateWin32Viewer : a class derived from G4OpenGLWin32Viewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLImmediateWin32Viewer.hh"

#include <GL/gl.h>

#include "G4ios.hh"

G4OpenGLImmediateWin32Viewer::G4OpenGLImmediateWin32Viewer
(G4OpenGLImmediateSceneHandler& scene,
 const G4String&  name):
G4OpenGLViewer (scene),
G4OpenGLWin32Viewer (scene),
G4OpenGLImmediateViewer (scene),
G4VViewer (scene, scene.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.
}

void G4OpenGLImmediateWin32Viewer::Initialise () {

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

void G4OpenGLImmediateWin32Viewer::DrawView () {

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
