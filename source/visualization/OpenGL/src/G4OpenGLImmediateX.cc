//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenGLImmediateX.cc,v 1.7 2001-07-11 10:08:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  27th March 1996
// OpenGL graphics system factory.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLImmediateXViewer.hh"

G4OpenGLImmediateX::G4OpenGLImmediateX ():
  G4VGraphicsSystem ("OpenGLImmediateX",
		     "OGLIX",
		     G4VisFeaturesOfOpenGLIX (),
		     G4VGraphicsSystem::threeD) {}

G4OpenGLImmediateX::~G4OpenGLImmediateX () {}

G4VSceneHandler* G4OpenGLImmediateX::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLImmediateSceneHandler (*this, name);
  G4cout << G4OpenGLImmediateSceneHandler::GetSceneCount ()
       << ' ' << fName << " scene handlers extanct." << G4endl;
  return    pScene;
}

G4VViewer* G4OpenGLImmediateX::CreateViewer (G4VSceneHandler& scene,
					 const G4String& name) {
  G4VViewer* pView =
    new G4OpenGLImmediateXViewer ((G4OpenGLImmediateSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateX::CreateViewer: error flagged by negative"
	" view id in G4OpenGLImmediateXViewer creation."
	"\n Destroying view and returning null pointer."
	   << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateX::CreateViewer: null pointer on"
      " new G4OpenGLImmediateXViewer." << G4endl;
  }
  return pView;
}

#endif
