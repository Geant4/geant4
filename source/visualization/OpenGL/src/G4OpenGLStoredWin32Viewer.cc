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
// $Id: G4OpenGLStoredWin32Viewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// 
// Class G4OpenGLStoredWin32Viewer : a class derived from G4OpenGLWin32Viewer and
//                             G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLStoredWin32Viewer.hh"

#include "G4OpenGLStoredSceneHandler.hh"

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

  glViewport (0, 0, getWinWidth(), getWinHeight());

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
    FinishView ();

  } else {

    // If kernel visit was needed, drawing and FinishView will already
    // have been done, so...
    if (!kernelVisitWasNeeded) {
      DrawDisplayLists ();
      FinishView ();
    } else {
    // However, union cutaways are implemented in DrawDisplayLists, so make
    // an extra pass...
      if (fVP.IsCutaway() &&
	  fVP.GetCutawayMode() == G4ViewParameters::cutawayUnion) {
	ClearView();
	DrawDisplayLists ();
	FinishView ();
      } else { // ADD TO AVOID KernelVisit=1 and nothing to display
        DrawDisplayLists ();
        FinishView ();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLStoredWin32Viewer::FinishView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fHDC) return;

  glFlush ();
  GLint renderMode;
  glGetIntegerv(GL_RENDER_MODE, &renderMode);
  if (renderMode == GL_RENDER) ::SwapBuffers(fHDC);

  // Empty the Windows message queue :
  MSG event;
  while ( ::PeekMessage(&event, NULL, 0, 0, PM_REMOVE) ) {
    ::TranslateMessage(&event);
    ::DispatchMessage (&event);
  }
}

#endif
