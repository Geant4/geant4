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
// $Id: G4XXX.cc,v 1.1 2001-08-17 23:04:33 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#include "G4XXX.hh"
#include "G4XXXSceneHandler.hh"
#include "G4XXXViewer.hh"

G4XXX::G4XXX():
  G4VGraphicsSystem("G4XXX",
		    "XXX",
		    "A template graphics driver",
		    G4VGraphicsSystem::noFunctionality) {}

G4XXX::~G4XXX() {}

G4VSceneHandler* G4XXX::CreateSceneHandler(const G4String& name) {
  G4VSceneHandler* pScene = new G4XXXSceneHandler(*this, name);
  G4cout << G4XXXSceneHandler::GetSceneCount()
         << ' ' << fName << " scene handlers extanct." << G4endl;
  return pScene;
}

G4VViewer* G4XXX::CreateViewer(G4VSceneHandler& scene,
			       const G4String& name) {
  G4VViewer* pView =
    new G4XXXViewer((G4XXXSceneHandler&) scene, name);
  if (pView) {
    if (pView->GetViewId() < 0) {
      G4cout <<
	"G4XXX::CreateViewer: ERROR flagged by negative"
        " view id in G4XXXViewer creation."
        "\n Destroying view and returning null pointer."
	     << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout <<
      "G4XXX::CreateViewer: ERROR: null pointer on new G4XXXViewer."
	   << G4endl;
  }
  return pView;
}
