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
// $Id: G4OpenGLImmediateWin32.cc,v 1.7 2002-02-24 01:47:58 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// OpenGLImmediateWin32 graphics system factory.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

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
       << ' ' << fName << " scene handlers extanct." << G4endl;
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
	   << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateWin32::CreateViewer: null pointer on"
      " new G4OpenGLImmediateWin32Viewer." << G4endl;
  }
  return pView;
}

#endif
