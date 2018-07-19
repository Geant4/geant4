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
// $Id: G4OpenGLImmediateXmViewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// 
// Andrew Walkden  10th February 1997
// Class G4OpenGLImmediateXmViewer : a class derived from G4OpenGLXmViewer
//                                     and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4OpenGLImmediateXmViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4ios.hh"
#include "G4Threading.hh"

G4OpenGLImmediateXmViewer::
G4OpenGLImmediateXmViewer(G4OpenGLImmediateSceneHandler& sceneHandler,
                          const G4String& name)
 : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
   G4OpenGLViewer (sceneHandler),
   G4OpenGLXmViewer (sceneHandler),
   G4OpenGLImmediateViewer (sceneHandler)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

// ensure a suitable window was found
  if (!vi_immediate) {
    G4cerr << "G4OpenGLImmediateXmViewer::G4OpenGLImmediateXmViewer -"
      " G4OpenGLXmViewer couldn't get a visual." << G4endl;  
    fViewId = -1;  // This flags an error.
    return;
  }
}

G4OpenGLImmediateXmViewer::~G4OpenGLImmediateXmViewer () {}

void G4OpenGLImmediateXmViewer::Initialise () {

  CreateGLXContext (vi_immediate);
  CreateMainWindow ();
  CreateFontLists ();

  InitializeGLView ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  glDrawBuffer (GL_FRONT);
}

void G4OpenGLImmediateXmViewer::DrawView () {

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateXmViewer::DrawView : \n");
#endif

  if(style!=G4ViewParameters::hlr &&
     haloing_enabled) {

    HaloingFirstPass ();
    NeedKernelVisit ();
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateXmViewer::DrawView : change param\n");
#endif
    ProcessView ();
    glFlush ();

    HaloingSecondPass ();

  }

#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateXmViewer::DrawView : need Kernel/Process/Finish\n");
#endif
  NeedKernelVisit ();  // Always need to visit G4 kernel.
  ProcessView ();
  FinishView ();

}

void G4OpenGLImmediateXmViewer::FinishView () {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateXmViewer::FinishView : \n");
#endif
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propogated before progressing.
  glFlush ();
}

#endif
