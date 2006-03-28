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
// $Id: G4XXXSG.cc,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  10th March 2005
// A template for a sophisticated graphics driver with a scene graph.
//?? Lines or sections marked like this require specialisation for your driver.

#ifdef G4VIS_BUILD_XXXSG_DRIVER

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
      G4cout <<
	"G4XXXSG::CreateViewer: ERROR flagged by negative"
        " view id in G4XXXSGViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout <<
      "G4XXXSG::CreateViewer: ERROR: null pointer on new G4XXXSGViewer."
	   << G4endl;
  }
  return pView;
}

#endif
