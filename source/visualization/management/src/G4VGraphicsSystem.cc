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
// $Id: G4VGraphicsSystem.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  27th March 1996
// Abstract interface class for graphics systems.

#include "G4VGraphicsSystem.hh"

#include "G4VisManager.hh"

G4VGraphicsSystem::~G4VGraphicsSystem () {}

G4VGraphicsSystem::G4VGraphicsSystem (const G4String& name,
				      Functionality f):
  fName (name),
  fNickname (""),
  fDescription (""),
  fFunctionality (f) {}

G4VGraphicsSystem::G4VGraphicsSystem (const G4String& name,
				      const G4String& nickname,
				      Functionality f):
  fName (name),
  fNickname (nickname),
  fDescription (""),
  fFunctionality (f) {}

G4VGraphicsSystem::G4VGraphicsSystem (const G4String& name,
				      const G4String& nickname,
				      const G4String& description,
				      Functionality f):
  fName (name),
  fNickname (nickname),
  fDescription (description),
  fFunctionality (f) {}

G4bool G4VGraphicsSystem::IsUISessionCompatible () const
{
  return true;
}

std::ostream& operator << (std::ostream& os, const G4VGraphicsSystem& gs) {
  G4VisManager* pVMan = G4VisManager::GetInstance ();
  const G4SceneHandlerList& scenes = pVMan -> GetAvailableSceneHandlers ();
  os << "Graphics System: " << gs.GetName ();
  if (gs.GetNickname () != "") {
    os << ", nickname: " << gs.GetNickname ();
  }
  if (gs.GetDescription () != "") {
    os << "\n  Description: " << gs.GetDescription ();
  }
  os << "\n  Functionality: " << G4int(gs.GetFunctionality());
  if (pVMan -> GetVerbosity() >= G4VisManager::parameters) {
    size_t nScenes = scenes.size ();
    if (nScenes) {
      G4int nScenesOfThisSystem = 0;
      for (size_t i = 0; i < nScenes; i++) {
	if (scenes [i] -> GetGraphicsSystem () == &gs) {
	  nScenesOfThisSystem++;
	}
      }
      if (nScenesOfThisSystem) {
	os << "\n  Its scenes are: ";
	for (size_t i = 0; i < nScenes; i++) {
	  if (scenes [i] -> GetGraphicsSystem () == &gs) {
	    os << "\n  " << *(scenes [i]);
	  }
	}
      }
      else {
	os << "\n  It has no scenes at present.";
      }
    }
    else {
      os << "\n  There are no scenes instantiated at present.";
    }
  }
  return os;
}
