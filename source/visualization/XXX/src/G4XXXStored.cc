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
// $Id: G4XXXStored.cc,v 1.1 2006-03-28 17:16:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  7th March 2006
// A template for a graphics driver with a store/database.
//?? Lines beginning like this require specialisation for your driver.

#include "G4XXXStored.hh"
#include "G4XXXStoredSceneHandler.hh"
#include "G4XXXStoredViewer.hh"

G4XXXStored::G4XXXStored():
  G4VGraphicsSystem("G4XXXStored",
		    "XXXStored",
		    "Graphics driver with a store/database",
		    G4VGraphicsSystem::threeD  //?? Your functionality
		    )
{}

G4XXXStored::~G4XXXStored() {}

G4VSceneHandler* G4XXXStored::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4XXXStoredSceneHandler(*this, name);
  return pScene;
}

G4VViewer* G4XXXStored::CreateViewer(G4VSceneHandler& scene,
			       const G4String& name) {
  G4VViewer* pView =
    new G4XXXStoredViewer((G4XXXStoredSceneHandler&) scene, name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cout <<
	"G4XXXStored::CreateViewer: ERROR flagged by negative"
        " view id in G4XXXStoredViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout <<
      "G4XXXStored::CreateViewer: ERROR: null pointer on new G4XXXStoredViewer."
	   << G4endl;
  }
  return pView;
}
