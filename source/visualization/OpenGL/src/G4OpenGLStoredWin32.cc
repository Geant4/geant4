// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredWin32.cc,v 1.1 1999-01-07 16:14:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// OpenGLStoredWin32 graphics system factory.


#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VScene.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLStoredWin32.hh"
#include "G4OpenGLStoredWin32View.hh"

G4OpenGLStoredWin32::G4OpenGLStoredWin32 ():
  G4VGraphicsSystem ("OpenGLStoredWin32",
		     "OGLSWin32",
		     G4VisFeaturesOfOpenGLSWin32 (),
		     G4VGraphicsSystem::threeD) {}

G4VScene* G4OpenGLStoredWin32::CreateScene () {
  G4VScene* pScene = new G4OpenGLStoredScene (*this);
  G4cout << G4OpenGLStoredScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VView* G4OpenGLStoredWin32::CreateView (G4VScene& scene) {
  G4VView* pView =
    new G4OpenGLStoredWin32View ((G4OpenGLStoredScene&) scene);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      delete pView;
      pView = 0;
      G4cerr << "G4OpenGLStoredWin32::CreateView: error flagged by"
	" negative view id in G4OpenGLStoredWin32View creation."
	"\n Destroying view and returning null pointer." << endl;
    }
  }
  else {
    G4cerr << "G4OpenGLStoredWin32::CreateView: null pointer on"
      " new G4OpenGLStoredWin32View." << endl;
  }
  return pView;
}

#endif

