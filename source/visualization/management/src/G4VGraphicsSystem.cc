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
// $Id: G4VGraphicsSystem.cc 99418 2016-09-21 09:18:42Z gcosmo $
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
  fDescription ("No description"),
  fFunctionality (f)
{
  fNicknames.push_back("No Nickname");
}

G4VGraphicsSystem::G4VGraphicsSystem (const G4String& name,
				      const G4String& nickname,
				      Functionality f):
  fName (name),
  fDescription ("No description"),
  fFunctionality (f)
{
  fNicknames.push_back(nickname);
}

G4VGraphicsSystem::G4VGraphicsSystem (const G4String& name,
				      const G4String& nickname,
				      const G4String& description,
				      Functionality f):
  fName (name),
  fDescription (description),
  fFunctionality (f)
{
  fNicknames.push_back(nickname);
}

G4bool G4VGraphicsSystem::IsUISessionCompatible () const
{
  return true;
}

std::ostream& operator << (std::ostream& os, const G4VGraphicsSystem& gs)
{
  os << "Graphics System: " << gs.GetName ();
  os << ", nicknames:"; for (const auto& nickname: gs.GetNicknames())
  {os << ' ' << nickname;}
  os << "\n  Description: " << gs.GetDescription ();
  os << "\n  Functionality: ";
  switch (gs.GetFunctionality()) {
    case G4VGraphicsSystem::noFunctionality:
      os << "None";
      break;
    case G4VGraphicsSystem::nonEuclidian:
      os << "nonEuclidian, e.g., tree representation of geometry hierarchy.";
      break;
    case G4VGraphicsSystem::twoD:
      os << "twoD: Simple 2D, e.g., X (no stored structures).";
      break;
    case G4VGraphicsSystem::twoDStore:
      os << "twoDStore: 2D with stored structures.";
      break;
    case G4VGraphicsSystem::threeD:
      os << "threeD: Passive 3D (with stored structures)";
      break;
    case G4VGraphicsSystem::threeDInteractive:
      os << "threeDInteractive: 3D with \"pick\" functionality.";
      break;
    case G4VGraphicsSystem::virtualReality:
      os << "virtualReality";
      break;
    case G4VGraphicsSystem::fileWriter:
      os << "fileWriter";
      break;
    default:
      os << "unknown";
      break;
  }

  G4VisManager* pVMan = G4VisManager::GetInstance ();
  const G4SceneHandlerList& scenes = pVMan -> GetAvailableSceneHandlers ();
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
