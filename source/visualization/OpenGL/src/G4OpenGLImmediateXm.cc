// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXm.cc,v 1.1 1999-01-07 16:14:57 gunter Exp $
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
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLImmediateXmView.hh"

G4OpenGLImmediateXm::G4OpenGLImmediateXm ():
  G4VGraphicsSystem ("OpenGLImmediateXm",
		     "OGLIXm",
		     G4VisFeaturesOfOpenGLIXm (),
		     G4VGraphicsSystem::threeD) {}

G4VScene* G4OpenGLImmediateXm::CreateScene (const G4String& name) {
  G4VScene* pScene = new G4OpenGLImmediateScene (*this, name);
  G4cout << G4OpenGLImmediateScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VView* G4OpenGLImmediateXm::CreateView (G4VScene& scene,
					  const G4String& name) {
  G4VView* pView =
    new G4OpenGLImmediateXmView ((G4OpenGLImmediateScene&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateXm::CreateView: error flagged by"
	" negative view id in G4OpenGLImmediateXmView creation."
	"\n Destroying view and returning null pointer." << endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateXm::CreateView: null pointer on"
      " new G4OpenGLImmediateXmView." << endl;
  }
  return pView;
}

#endif
