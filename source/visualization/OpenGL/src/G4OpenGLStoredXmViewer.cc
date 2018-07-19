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
// $Id: G4OpenGLStoredXmViewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// 
// Andrew Walkden  10th February 1997
// Class G4OpenGLStoredXmViewer : a class derived from G4OpenGLXmViewer 
//                              and G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLStoredXmViewer.hh"

#include "G4OpenGLStoredSceneHandler.hh"
#include "G4ios.hh"
#include "G4Threading.hh"

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

G4OpenGLStoredXmViewer::~G4OpenGLStoredXmViewer () {
  GetSceneHandler()->RemoveViewerFromList(this);
}

void G4OpenGLStoredXmViewer::Initialise () {

  CreateGLXContext (vi_stored);
  CreateMainWindow ();
  CreateFontLists ();

  InitializeGLView ();

  glDrawBuffer (GL_BACK);
}

void G4OpenGLStoredXmViewer::DrawView () {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLStoredXmViewer::DrawView \n");
#endif

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
#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLStoredXmViewer::DrawView () flush\n");
#endif
    glFlush ();

    HaloingSecondPass ();

    DrawDisplayLists ();
    FinishView ();

  } else {

#ifdef G4DEBUG_VIS_OGL
    printf("G4OpenGLStoredXmViewer::DrawView not hlr \n");
#endif
    // If kernel visit was needed, drawing and FinishView will already
    // have been done, so...
    if (!kernelVisitWasNeeded) {
#ifdef G4DEBUG_VIS_OGL
      printf("G4OpenGLStoredXmViewer::ComputeView Don't need kernel Visit \n");
#endif
      DrawDisplayLists ();
      FinishView ();
    } else {
#ifdef G4DEBUG_VIS_OGL
      printf("G4OpenGLStoredXmViewer::ComputeView Need kernel Visit \n");
#endif
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
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLStoredXmViewer::DrawView  ^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
#endif
}

void G4OpenGLStoredXmViewer::FinishView () {
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propogated before progressing.

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLStoredXmViewer::FinishView () flush \n");
#endif
  glFlush (); //FIXME

  GLint renderMode;
  glGetIntegerv(GL_RENDER_MODE, &renderMode);
  if (renderMode == GL_RENDER) glXSwapBuffers (dpy, win);  
}

#endif
