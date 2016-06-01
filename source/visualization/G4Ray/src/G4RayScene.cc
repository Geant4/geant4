// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayScene.cc,v 2.3 1998/11/06 13:41:55 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing scene.

#include "G4RayScene.hh"

G4RayScene::G4RayScene (G4VGraphicsSystem& system,
			const G4String& name):
G4VScene (system, fSceneIdCount++, name)
{
  fSceneCount++;
}

G4RayScene::~G4RayScene ()
{
  fSceneCount--;
}

G4int G4RayScene::GetSceneCount () {
  return fSceneCount;
}

void G4RayScene::AddPrimitive (const G4Polyline&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;
    G4cout << "  void G4RayScene::AddPrimitive (const G4Polyline&) called."
	   << "  \nG4Ray cannot handle polylines."
	   << endl;
  }
}

void G4RayScene::AddPrimitive (const G4Text&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;
    G4cout << "  void G4RayScene::AddPrimitive (const G4Text&) called."
	   << "  \nG4Ray cannot handle text."
	   << endl;
  }
}
  
void G4RayScene::AddPrimitive (const G4Circle&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;
    G4cout << "  void G4RayScene::AddPrimitive (const G4Circle&) called."
	   << "  \nG4Ray cannot handle markers."
	   << endl;
  }
}

void G4RayScene::AddPrimitive (const G4Square&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;
    G4cout << "  void G4RayScene::AddPrimitive (const G4Square&) called."
	   << "  \nG4Ray cannot handle markers."
	   << endl;
  }
}

void G4RayScene::AddPrimitive (const G4Polyhedron&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;  
    G4cout << "void G4RayScene::AddPrimitive (const G4Polyhedron&) called."
	   << "  \nG4Ray cannot handle polyhedra."
	   << endl;
  }
}

void G4RayScene::AddPrimitive (const G4NURBS&) {
  static G4bool warned = false;
  if (!warned) {
    warned = true;
    G4cout << "***WARNING***" << endl;
    G4cout << "void G4RayScene::AddPrimitive (const G4NURBS&) called."
	   << "  \nG4Ray cannot handle NURBS."
	 << endl;
  }
}

G4int G4RayScene::fSceneIdCount = 0;

G4int G4RayScene::fSceneCount = 0;
