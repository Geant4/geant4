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
// $Id: G4OpenGLStoredXm.cc 91686 2015-07-31 09:40:08Z gcosmo $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLStoredSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLStoredXm.hh"
#include "G4OpenGLStoredXmViewer.hh"
#include "G4OpenGLViewerMessenger.hh"
#include "G4OpenGLXmViewerMessenger.hh"

G4OpenGLStoredXm::G4OpenGLStoredXm ():
  G4OpenGLXm ("OpenGLStoredXm",
              "OGLSXm",
              G4VisFeaturesOfOpenGLSXm (),
              G4VGraphicsSystem::threeD)
{
  G4OpenGLViewerMessenger::GetInstance();
  G4OpenGLXmViewerMessenger::GetInstance();
}

G4VSceneHandler* G4OpenGLStoredXm::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLStoredSceneHandler (*this, name);
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
