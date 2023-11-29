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
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#include "G4ASCIITree.hh"
#include "G4ASCIITreeSceneHandler.hh"
#include "G4ASCIITreeViewer.hh"
#include "G4ASCIITreeMessenger.hh"

#define G4warn G4cout

G4ASCIITree::G4ASCIITree ():
  G4VTree ("ASCIITree",
	   "ATree",
	   "A graphics system to dump geometry hierarchy"
	   "\n  to standard output as an ASCII stream.",
	   G4VGraphicsSystem::nonEuclidian),
  fVerbosity(1),
  fOutFileName ("G4cout")
{
  fpMessenger = new G4ASCIITreeMessenger(this);
}

G4ASCIITree::~G4ASCIITree () {
  delete fpMessenger;
}

G4VSceneHandler* G4ASCIITree::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4ASCIITreeSceneHandler (*this, name);
  return pScene;
}

G4VViewer* G4ASCIITree::CreateViewer (G4VSceneHandler& scene,
				  const G4String& name) {
  G4VViewer* pView =
    new G4ASCIITreeViewer ((G4ASCIITreeSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4warn << "G4ASCIITree::CreateViewer: ERROR flagged by negative"
        " view id in G4ASCIITreeViewer creation."
        "\n Destroying view and returning null pointer."
           << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4warn << "G4ASCIITree::CreateViewer: ERROR: null pointer on"
      " new G4ASCIITreeViewer." << G4endl;
  }
  return pView;
}
