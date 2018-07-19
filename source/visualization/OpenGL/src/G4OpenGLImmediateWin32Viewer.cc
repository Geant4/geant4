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
// $Id: G4OpenGLImmediateWin32Viewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// 
// Class G4OpenGLImmediateWin32Viewer : a class derived from G4OpenGLWin32Viewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include "G4OpenGLImmediateWin32Viewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4ios.hh"

G4OpenGLImmediateWin32Viewer::G4OpenGLImmediateWin32Viewer
(G4OpenGLImmediateSceneHandler& sceneHandler,
 const G4String&  name):
G4OpenGLViewer (sceneHandler),
G4OpenGLWin32Viewer (sceneHandler),
G4OpenGLImmediateViewer (sceneHandler),
G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name) {

  if (fViewId < 0) return;  // In case error in base class instantiation.
}

void G4OpenGLImmediateWin32Viewer::Initialise () {

// ensure a suitable window was found

  CreateGLWin32Context ();
  CreateMainWindow ();
  CreateFontLists ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  glDrawBuffer (GL_FRONT);

  // clear the buffers and window.
  ClearView ();
  FinishView ();

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
}

void G4OpenGLImmediateWin32Viewer::DrawView () {

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  glDrawBuffer (GL_FRONT);

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

  glViewport (0, 0, getWinWidth(), getWinHeight());

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

//////////////////////////////////////////////////////////////////////////////
void G4OpenGLImmediateWin32Viewer::FinishView (
) 
//////////////////////////////////////////////////////////////////////////////
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
{
  if(!fHDC) return;

  glFlush ();

  // Empty the Windows message queue :
  MSG event;
  while ( ::PeekMessage(&event, NULL, 0, 0, PM_REMOVE) ) {
    ::TranslateMessage(&event);
    ::DispatchMessage (&event);
  }
}

#endif
