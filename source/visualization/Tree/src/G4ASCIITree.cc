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
// $Id: G4ASCIITree.cc,v 1.7 2001-08-05 19:02:12 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy.

#include "G4ASCIITree.hh"
#include "G4ASCIITreeSceneHandler.hh"
#include "G4ASCIITreeViewer.hh"
#include "G4ASCIITreeMessenger.hh"

G4ASCIITree::G4ASCIITree ():
  G4VTree ("ASCIITree",
	   "ATree",
	   "A graphics system to dump geometry hierarchy"
	   "\n  to standard output as an ASCII stream.",
	   G4VGraphicsSystem::nonEuclidian),
  fVerbosity(0)
{
  fpMessenger = new G4ASCIITreeMessenger(this);
}

G4ASCIITree::~G4ASCIITree () {
  delete fpMessenger;
}

G4VSceneHandler* G4ASCIITree::CreateSceneHandler (const G4String& name) {
  G4VSceneHandler* pScene = new G4ASCIITreeSceneHandler (*this, name);
  G4cout << G4ASCIITreeSceneHandler::GetSceneCount ()
	 << ' ' << fName << " scene handlers extanct." << G4endl;
  return pScene;
}

G4VViewer* G4ASCIITree::CreateViewer (G4VSceneHandler& scene,
				  const G4String& name) {
  G4VViewer* pView =
    new G4ASCIITreeViewer ((G4ASCIITreeSceneHandler&) scene, name);
  if (pView) {
    if (pView -> GetViewId () < 0) {
      G4cout << "G4ASCIITree::CreateViewer: ERROR flagged by negative"
        " view id in G4ASCIITreeViewer creation."
        "\n Destroying view and returning null pointer."
           << G4endl;
      delete pView;
      pView = 0;
    }
  }
  else {
    G4cout << "G4ASCIITree::CreateViewer: ERROR: null pointer on"
      " new G4ASCIITreeViewer." << G4endl;
  }
  return pView;
}
