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
// $Id: G4OpenGLStoredXm.cc,v 1.7 2002-02-24 01:48:20 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredXm.hh"
#include "G4OpenGLStoredXmViewer.hh"

G4OpenGLStoredXm::G4OpenGLStoredXm ():
  G4VGraphicsSystem ("OpenGLStoredXm",
		     "OGLSXm",
		     G4VisFeaturesOfOpenGLSXm (),
		     G4VGraphicsSystem::threeD) {}

G4VSceneHandler* G4OpenGLStoredXm::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLStoredSceneHandler (*this, name);
  G4cout << G4OpenGLStoredSceneHandler::GetSceneCount ()
       << ' ' << fName << " scene handlers extanct." << G4endl;
  return    pScene;
}

G4VViewer* G4OpenGLStoredXm::CreateViewer (G4VSceneHandler& scene, const G4String& name) {
  G4VViewer* pView =
    new G4OpenGLStoredXmViewer ((G4OpenGLStoredSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      delete pView;
      pView = 0;
      G4cerr << "G4OpenGLStoredXm::CreateViewer: error flagged by"
	" negative view id in G4OpenGLStoredXmViewer creation."
	"\n Destroying view and returning null pointer." << G4endl;
    }
  }
  else {
    G4cerr << "G4OpenGLStoredXm::CreateViewer: null pointer on"
      " new G4OpenGLStoredXmViewer." << G4endl;
  }
  return pView;
}

#endif
