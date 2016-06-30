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
// $Id: G4OpenGLImmediateXViewer.cc 97241 2016-05-30 12:06:54Z gcosmo $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateXViewer : a class derived from G4OpenGLXViewer and
//                                G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include "G4OpenGLImmediateXViewer.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4ios.hh"
#include "G4Threading.hh"

G4OpenGLImmediateXViewer::
G4OpenGLImmediateXViewer (G4OpenGLImmediateSceneHandler& sceneHandler,
			  const G4String&  name)
 : G4VViewer (sceneHandler, sceneHandler.IncrementViewCount (), name),
   G4OpenGLViewer (sceneHandler),
   G4OpenGLXViewer (sceneHandler),
   G4OpenGLImmediateViewer (sceneHandler)
{
  if (fViewId < 0) return;  // In case error in base class instantiation.

// ensure a suitable window was found
  if (!vi_immediate) {
    G4cerr << "G4OpenGLImmediateXViewer::G4OpenGLImmediateXViewer -"
      " G4OpenGLXViewer couldn't get a visual." << G4endl;  
    fViewId = -1;  // This flags an error.
    return;
  }
}

G4OpenGLImmediateXViewer::~G4OpenGLImmediateXViewer () {}

void G4OpenGLImmediateXViewer::Initialise () {

  CreateGLXContext (vi_immediate);
  CreateMainWindow ();
  CreateFontLists ();

  InitializeGLView ();

  // If a double buffer context has been forced upon us, ignore the
  // back buffer for this OpenGLImmediate view.
  glDrawBuffer (GL_FRONT);

  glDepthFunc (GL_LEQUAL);
  glDepthMask (GL_TRUE);
}

void G4OpenGLImmediateXViewer::DrawView () {

  G4ViewParameters::DrawingStyle style = GetViewParameters().GetDrawingStyle();

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

void G4OpenGLImmediateXViewer::FinishView () {
  glXWaitGL (); //Wait for effects of all previous OpenGL commands to
                //be propogated before progressing.
  glFlush ();
}

#endif
