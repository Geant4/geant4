//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenGLImmediateWin32.cc 66373 2012-12-18 09:41:34Z gcosmo $
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
#include "G4OpenGLViewerMessenger.hh"
#include "G4OpenGLImmediateSceneHandler.hh"

G4OpenGLImmediateWin32::G4OpenGLImmediateWin32 ():
  G4VGraphicsSystem ("OpenGLImmediateWin32",
		     "OGLIWin32",
		     G4VisFeaturesOfOpenGLIWin32 (),
		     G4VGraphicsSystem::threeD)
{
  G4OpenGLViewerMessenger::GetInstance();
}

G4VSceneHandler* G4OpenGLImmediateWin32::CreateSceneHandler
(const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLImmediateSceneHandler (*this, name);
  return    pScene;
}

G4VViewer* G4OpenGLImmediateWin32::CreateViewer
(G4VSceneHandler& scene, const G4String& name) {
  G4VViewer* pView =
    new G4OpenGLImmediateWin32Viewer
    ((G4OpenGLImmediateSceneHandler&) scene, name);
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
