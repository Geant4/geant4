// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredX.cc,v 1.6 2000-08-19 18:34:36 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.


#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredX.hh"
#include "G4OpenGLStoredXViewer.hh"

G4OpenGLStoredX::G4OpenGLStoredX ():
  G4VGraphicsSystem ("OpenGLStoredX",
		     "OGLSX",
		     G4VisFeaturesOfOpenGLSX (),
		     G4VGraphicsSystem::threeD) {}

G4OpenGLStoredX::~G4OpenGLStoredX () {}

G4VSceneHandler* G4OpenGLStoredX::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLStoredSceneHandler (*this, name);
  G4cout << G4OpenGLStoredSceneHandler::GetSceneCount ()
       << ' ' << fName << " scene handlers extanct." << G4endl;
  return    pScene;
}

G4VViewer* G4OpenGLStoredX::CreateViewer (G4VSceneHandler& scene,
				      const G4String& name) {
  G4VViewer* pView =
    new G4OpenGLStoredXViewer ((G4OpenGLStoredSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      delete pView;
      pView = 0;
      G4cerr << "G4OpenGLStoredX::CreateViewer: error flagged by"
	" negative view id in G4OpenGLStoredXViewer creation."
	"\n Destroying view and returning null pointer." << G4endl;
    }
  }
  else {
    G4cerr << "G4OpenGLStoredX::CreateViewer: null pointer on"
      " new G4OpenGLStoredXViewer." << G4endl;
  }
  return pView;
}

#endif

