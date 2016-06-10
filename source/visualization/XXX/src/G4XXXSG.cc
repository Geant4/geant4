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
// $Id: G4XXXSG.cc 85582 2014-10-31 09:07:30Z gcosmo $
//
// 
// John Allison  10th March 2005
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines or sections marked like this require specialisation for your driver.

#include "G4XXXSG.hh"
#include "G4XXXSGSceneHandler.hh"
#include "G4XXXSGViewer.hh"

G4XXXSG::G4XXXSG():
  G4VGraphicsSystem("G4XXXSG",
		    "XXXSG",
		    "Graphics driver with scene graph",
		    G4VGraphicsSystem::threeD) {}  //?? Your functionality.

G4XXXSG::~G4XXXSG() {}

G4VSceneHandler* G4XXXSG::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4XXXSGSceneHandler(*this, name);
  return pScene;
}

G4VViewer* G4XXXSG::CreateViewer(G4VSceneHandler& scene,
			       const G4String& name) {
  G4VViewer* pView =
    new G4XXXSGViewer((G4XXXSGSceneHandler&) scene, name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cerr <<
	"G4XXXSG::CreateViewer: ERROR flagged by negative"
        " view id in G4XXXSGViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr <<
      "G4XXXSG::CreateViewer: ERROR: null pointer on new G4XXXSGViewer."
	   << G4endl;
  }
  return pView;
}
