// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateWin32.cc,v 1.3 1999-01-11 00:47:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// OpenGLImmediateWin32 graphics system factory.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLImmediateWin32Viewer.hh"

G4OpenGLImmediateWin32::G4OpenGLImmediateWin32 ():
  G4VGraphicsSystem ("OpenGLImmediateWin32",
		     "OGLIWin32",
		     G4VisFeaturesOfOpenGLIWin32 (),
		     G4VGraphicsSystem::threeD) {}

G4VSceneHandler* G4OpenGLImmediateWin32::CreateSceneHandler () {
  G4VSceneHandler* pScene = new G4OpenGLImmediateSceneHandler (*this);
  G4cout << G4OpenGLImmediateSceneHandler::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VViewer* G4OpenGLImmediateWin32::CreateViewer (G4VSceneHandler& scene) {
  G4VViewer* pView =
    new G4OpenGLImmediateWin32Viewer ((G4OpenGLImmediateSceneHandler&) scene);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateWin32::CreateViewer: error flagged by negative"
	" view id in G4OpenGLImmediateWin32Viewer creation."
	"\n Destroying view and returning null pointer."
	   << endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateWin32::CreateViewer: null pointer on"
      " new G4OpenGLImmediateWin32Viewer." << endl;
  }
  return pView;
}

#endif
