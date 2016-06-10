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
// $Id: G4OpenGLImmediateWt.cc 75567 2013-11-04 11:35:11Z gcosmo $
//
// 
// OpenGLImmediateWt graphics system factory.

#ifdef G4VIS_BUILD_OPENGLWT_DRIVER

#include "G4VisFeaturesOfOpenGL.hh"
#include "G4VSceneHandler.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLImmediateWt.hh"
#include "G4OpenGLImmediateWtViewer.hh"
#include "G4OpenGLViewerMessenger.hh"
#include "G4OpenGLImmediateSceneHandler.hh"
#include "G4UIWt.hh"
#include "G4UImanager.hh"

#include <Wt/WLabel>
#include <Wt/WContainerWidget>

G4OpenGLImmediateWt::G4OpenGLImmediateWt ():
  G4VGraphicsSystem ("OpenGLImmediateWt",
		     "OGLIWt",
		     G4VisFeaturesOfOpenGLIWt (),
		     G4VGraphicsSystem::threeD)
{
  G4OpenGLViewerMessenger::GetInstance();
}

G4VSceneHandler* G4OpenGLImmediateWt::CreateSceneHandler
(const G4String& name) {
  G4VSceneHandler* pScene = new G4OpenGLImmediateSceneHandler (*this, name);
  return    pScene;
}

G4VViewer* G4OpenGLImmediateWt::CreateViewer
(G4VSceneHandler& scene, const G4String& name) {
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWt::CreateViewer \n");
#endif
  G4VViewer* pView = 0;
  // G4VisManager::CreateViewer
  // -> G4OpenGLImmediateWt::CreateViewer
  //   -> G4OpenGLImmediateWtViewer
  // !! ImmediateWtViewer should be the container
  
  // G4VisManager::CreateViewer
  // G4OpenGLImmediateQtViewer::Initialise()
  // ->G4OpenGLQtViewer::CreateMainWindow(container, name)
  //   fUiQt->AddTabWidget  (container fWindow)

  // Create a container for the widget
  fWGLContainer = new Wt::WContainerWidget();
  
  pView = new G4OpenGLImmediateWtViewer
  ((G4OpenGLImmediateSceneHandler&) scene, fWGLContainer, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4OpenGLImmediateWt::CreateViewer: error flagged by negative"
      " view id in G4OpenGLImmediateWtViewer creation."
      "\n Destroying view and returning null pointer."
      << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4OpenGLImmediateWt::CreateViewer: null pointer on"
    " new G4OpenGLImmediateWtViewer." << G4endl;
  }
#ifdef G4DEBUG_VIS_OGL
  printf("G4OpenGLImmediateWt::CreateViewer END \n");
#endif
  return pView;
}


#endif
