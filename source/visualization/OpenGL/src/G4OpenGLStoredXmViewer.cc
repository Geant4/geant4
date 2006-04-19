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
// $Id: G4OpenGLStoredXmViewer.cc,v 1.17 2006-04-19 12:02:37 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// Class G4OpenGLStoredXmViewer : a class derived from G4OpenGLXmViewer 
//                              and G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLStoredXmViewer.hh"

#include "G4ios.hh"

G4OpenGLStoredXmViewer::
G4OpenGLStoredXmViewer (G4OpenGLStoredSceneHandler& sceneHandler,
			const G4String& name)
 : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
   G4OpenGLViewer (sceneHandler),
   G4OpenGLXmViewer (sceneHandler),
   G4OpenGLStoredViewer (sceneHandler)
{

  if (fViewId < 0) return;  // In case error in base class instantiation.

  if (!vi_stored) {
    fViewId = -1;  // This flags an error.
    G4cerr << "G4OpenGLStoredXmViewer::G4OpenGLStoredXmViewer -"
      " G4OpenGLXmViewer couldn't get a visual." << G4endl;
    return;
  }
}

G4OpenGLStoredXmViewer::~G4OpenGLStoredXmViewer () {}

void G4OpenGLStoredXmViewer::Initialise () {

  CreateGLXContext (vi_stored);
  CreateMainWindow ();
  CreateFontLists ();

  InitializeGLView ();

// clear the buffers and window.
  ClearView ();
  FinishView ();
  
  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
  
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void G4OpenGLStoredXmViewer::DrawView () {

  //Make sure current viewer is attached and clean...
  glXMakeCurrent (dpy, win, cx);
  glViewport (0, 0, WinSize_x, WinSize_y);

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  //See if things have changed from last time and remake if necessary...
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
