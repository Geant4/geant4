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
// Satoshi Tanaka 31th May 2001
// A graphics system to dump geometry hierarchy to GAG.

#include "G4GAGTree.hh"
#include "G4GAGTreeSceneHandler.hh"
#include "G4GAGTreeViewer.hh"
#include "G4GAGTreeMessenger.hh"

G4GAGTree::G4GAGTree ():
  G4VTree ("GAGTree",
	   "GAGTree",
	   "A graphics system to dump geometry hierarchy"
	   "\n  to GAG.",
	   G4VGraphicsSystem::nonEuclidian) {
  fpMessenger = new G4GAGTreeMessenger(this);
}

G4GAGTree::~G4GAGTree () {
  delete fpMessenger;
}

G4VSceneHandler* G4GAGTree::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4GAGTreeSceneHandler (*this, name);
  G4cout << G4GAGTreeSceneHandler::GetSceneCount ()
	 << ' ' << fName << " scene handlers extanct." << G4endl;
  return pScene;
}

G4VViewer* G4GAGTree::CreateViewer (G4VSceneHandler& scene,
				  const G4String& name) {
  G4VViewer* pView =
    new G4GAGTreeViewer ((G4GAGTreeSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cerr << "G4GAGTree::CreateViewer: error flagged by negative"
        " view id in G4GAGTreeViewer creation."
        "\n Destroying view and returning null pointer."
           << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cerr << "G4GAGTree::CreateViewer: null pointer on"
      " new G4GAGTreeViewer." << G4endl;
  }
  return pView;
}
