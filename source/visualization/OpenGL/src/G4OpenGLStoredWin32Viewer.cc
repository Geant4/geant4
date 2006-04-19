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
// $Id: G4OpenGLStoredWin32Viewer.cc,v 1.14 2006-04-19 12:02:37 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLStoredWin32Viewer : a class derived from G4OpenGLWin32Viewer and
//                             G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLStoredWin32Viewer.hh"

#include "G4ios.hh"

G4OpenGLStoredWin32Viewer::G4OpenGLStoredWin32Viewer
(G4OpenGLStoredSceneHandler& sceneHandler,
 const G4String&  name):
G4OpenGLViewer (sceneHandler),
G4OpenGLWin32Viewer (sceneHandler),
G4OpenGLStoredViewer (sceneHandler),
G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.
}

void G4OpenGLStoredWin32Viewer::Initialise () {

  //Check that G4OpenGLWin32Viewer got a double buffered colour visual

  CreateGLWin32Context ();
  CreateMainWindow ();
  CreateFontLists ();

// clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
}

void G4OpenGLStoredWin32Viewer::DrawView () {

  //Make sure current viewer is attached and clean...
  //Win32 version needed
  //  glXMakeCurrent (dpy, win, cx);
  glViewport (0, 0, WinSize_x, WinSize_y);

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  //See if things have changed from last time and remake if necessary...
  // The fNeedKernelVisit flag might have been set by the user in
  // /vis/viewer/rebuild, but if not, make decision and set flag only
  // if necessary...
  if (!fNeedKernelVisit) KernelVisitDecision ();
  G4bool kernelVisitWasNeeded = fNeedKernelVisit; // Keep (ProcessView resets).
  ProcessView ();

  if(style!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    DrawDisplayLists ();
    glFlush ();

    HaloingSecondPass ();

    DrawDisplayLists ();

  } else {

    // If kernel visit was needed, it will already have been drawn, so...
    if (!kernelVisitWasNeeded) DrawDisplayLists ();

  }

  FinishView ();
  
}

#endif
