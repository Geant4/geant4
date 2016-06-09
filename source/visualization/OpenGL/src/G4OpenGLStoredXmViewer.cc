//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLStoredXmViewer.cc,v 1.19 2006/06/29 21:19:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
    FinishView ();

  } else {

    // If kernel visit was needed, drawing and FinishView will already
    // have been done, so...
    if (!kernelVisitWasNeeded) {
      DrawDisplayLists ();
      FinishView ();
    }
  }
}

#endif
