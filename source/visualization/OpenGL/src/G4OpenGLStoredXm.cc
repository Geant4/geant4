// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredXm.cc,v 1.1 1999-01-07 16:14:59 gunter Exp $
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
#include "G4VScene.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLStoredXm.hh"
#include "G4OpenGLStoredXmView.hh"

G4OpenGLStoredXm::G4OpenGLStoredXm ():
  G4VGraphicsSystem ("OpenGLStoredXm",
		     "OGLSXm",
		     G4VisFeaturesOfOpenGLSXm (),
		     G4VGraphicsSystem::threeD) {}

G4VScene* G4OpenGLStoredXm::CreateScene (const G4String& name) {
  G4VScene* pScene = new G4OpenGLStoredScene (*this, name);
  G4cout << G4OpenGLStoredScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VView* G4OpenGLStoredXm::CreateView (G4VScene& scene, const G4String& name) {
  G4VView* pView =
    new G4OpenGLStoredXmView ((G4OpenGLStoredScene&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      delete pView;
      pView = 0;
      G4cerr << "G4OpenGLStoredXm::CreateView: error flagged by"
	" negative view id in G4OpenGLStoredXmView creation."
	"\n Destroying view and returning null pointer." << endl;
    }
  }
  else {
    G4cerr << "G4OpenGLStoredXm::CreateView: null pointer on"
      " new G4OpenGLStoredXmView." << endl;
  }
  return pView;
}

#endif
