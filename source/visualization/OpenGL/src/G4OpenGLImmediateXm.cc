// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXm.cc,v 1.6 2001-02-03 18:39:29 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLImmediateXmViewer.hh"

G4OpenGLImmediateXm::G4OpenGLImmediateXm ():
  G4VGraphicsSystem ("OpenGLImmediateXm",
		     "OGLIXm",
		     G4VisFeaturesOfOpenGLIXm (),
		     G4VGraphicsSystem::threeD) {}

G4OpenGLImmediateXm::~G4OpenGLImmediateXm () {}

G4VSceneHandler* G4OpenGLImmediateXm::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLImmediateSceneHandler (*this, name);
  G4cout << G4OpenGLImmediateSceneHandler::GetSceneCount ()
       << ' ' << fName << " scene handlers extanct." << G4endl;
  return    pScene;
}

G4VViewer* G4OpenGLImmediateXm::CreateViewer (G4VSceneHandler& scene,
					  const G4String& name) {
  G4VViewer* pView =
    new G4OpenGLImmediateXmViewer ((G4OpenGLImmediateSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateXm::CreateViewer: error flagged by"
	" negative view id in G4OpenGLImmediateXmViewer creation."
	"\n Destroying view and returning null pointer." << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateXm::CreateViewer: null pointer on"
      " new G4OpenGLImmediateXmViewer." << G4endl;
  }
  return pView;
}

#endif
