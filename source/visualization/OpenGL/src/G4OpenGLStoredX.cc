// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredX.cc,v 1.1 1999-01-07 16:14:58 gunter Exp $
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
#include "G4VScene.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLStoredX.hh"
#include "G4OpenGLStoredXView.hh"

G4OpenGLStoredX::G4OpenGLStoredX ():
  G4VGraphicsSystem ("OpenGLStoredX",
		     "OGLSX",
		     G4VisFeaturesOfOpenGLSX (),
		     G4VGraphicsSystem::threeD) {}

G4VScene* G4OpenGLStoredX::CreateScene (const G4String& name) {
  G4VScene* pScene = new G4OpenGLStoredScene (*this, name);
  G4cout << G4OpenGLStoredScene::GetSceneCount ()
       << ' ' << fName << " scenes extanct." << endl;
  return    pScene;
}

G4VView* G4OpenGLStoredX::CreateView (G4VScene& scene,
				      const G4String& name) {
  G4VView* pView =
    new G4OpenGLStoredXView ((G4OpenGLStoredScene&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      delete pView;
      pView = 0;
      G4cerr << "G4OpenGLStoredX::CreateView: error flagged by"
	" negative view id in G4OpenGLStoredXView creation."
	"\n Destroying view and returning null pointer." << endl;
    }
  }
  else {
    G4cerr << "G4OpenGLStoredX::CreateView: null pointer on"
      " new G4OpenGLStoredXView." << endl;
  }
  return pView;
}

#endif

