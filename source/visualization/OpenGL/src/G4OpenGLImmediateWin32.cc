// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateWin32.cc,v 1.1 1999-01-07 16:14:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// OpenGLImmediateWin32 graphics system factory.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VScene.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLImmediateWin32.hh"
#include "G4OpenGLImmediateWin32View.hh"

G4OpenGLImmediateWin32::G4OpenGLImmediateWin32 ():
  G4VGraphicsSystem ("OpenGLImmediateWin32",
		     "OGLIWin32",
		     G4VisFeaturesOfOpenGLIWin32 (),
		     G4VGraphicsSystem::threeD) {}

G4VScene* G4OpenGLImmediateWin32::CreateScene () {
  G4VScene* pScene = new G4OpenGLImmediateScene (*this);
  G4cout << G4OpenGLImmediateScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VView* G4OpenGLImmediateWin32::CreateView (G4VScene& scene) {
  G4VView* pView =
    new G4OpenGLImmediateWin32View ((G4OpenGLImmediateScene&) scene);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateWin32::CreateView: error flagged by negative"
	" view id in G4OpenGLImmediateWin32View creation."
	"\n Destroying view and returning null pointer."
	   << endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateWin32::CreateView: null pointer on"
      " new G4OpenGLImmediateWin32View." << endl;
  }
  return pView;
}

#endif
