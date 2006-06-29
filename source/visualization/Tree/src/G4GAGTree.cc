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
  return pScene;
}

G4VViewer* G4GAGTree::CreateViewer (G4VSceneHandler& scene,
				  const G4String& name) {
  G4VViewer* pView =
    new G4GAGTreeViewer ((G4GAGTreeSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cout << "G4GAGTree::CreateViewer: ERROR flagged by negative"
        " view id in G4GAGTreeViewer creation."
        "\n Destroying view and returning null pointer."
           << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout << "G4GAGTree::CreateViewer: ERROR: null pointer on"
      " new G4GAGTreeViewer." << G4endl;
  }
  return pView;
}
