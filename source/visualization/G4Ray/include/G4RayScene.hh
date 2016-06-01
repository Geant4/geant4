// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayScene.hh,v 2.2 1998/11/06 13:41:50 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing scene.

#ifndef G4RAYSCENE_HH
#define G4RAYSCENE_HH

#include "G4VScene.hh"

class G4RayScene: public G4VScene {

public:
  G4RayScene (G4VGraphicsSystem& system, const G4String& name = "");
  ~G4RayScene ();
  static G4int GetSceneCount ();

  // The following functions don't make sense for G4Ray...
  void AddPrimitive (const G4Polyline&);
  void AddPrimitive (const G4Text&);
  void AddPrimitive (const G4Circle&);
  void AddPrimitive (const G4Square&);
  void AddPrimitive (const G4Polymarker&);
  void AddPrimitive (const G4Polyhedron&);
  void AddPrimitive (const G4NURBS&);
  void AddThis (const G4Box&);
  void AddThis (const G4Cons&);
  void AddThis (const G4Tubs&);
  void AddThis (const G4Trd&);
  void AddThis (const G4Trap&);
  void AddThis (const G4Sphere&);
  void AddThis (const G4Para&);
  void AddThis (const G4Torus&);
  void AddThis (const G4VSolid&);

private:
  static G4int    fSceneIdCount;  // static counter for Ray scenes.
  static G4int    fSceneCount;    // No. of extanct scenes.
};

#include "G4RayScene.icc"

#endif
