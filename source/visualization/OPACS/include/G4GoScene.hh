// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GoScene.hh,v 2.1 1998/11/06 13:42:01 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo scene - creates Wo Display lists.

#ifndef G4GOSCENE_HH
#define G4GOSCENE_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

//Go
#include <ONode.h>
#include <OColormap.h>
//G4
#include "globals.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VScene.hh"

class G4GoScene: public G4VScene {

public:
  G4GoScene            (G4VGraphicsSystem& system,
			const G4String& name = "");
  ~G4GoScene           ();
  void AddPrimitive    (const G4Polyline&);
  void AddPrimitive    (const G4Text&);
  void AddPrimitive    (const G4Circle&);
  void AddPrimitive    (const G4Square&);
  void AddPrimitive    (const G4Polymarker&);
  void AddPrimitive    (const G4Polyhedron&);
  void AddPrimitive    (const G4NURBS&);

  void AddThis         (const G4Box&);
  void AddThis         (const G4Cons&);
  void AddThis         (const G4Tubs&);
  void AddThis         (const G4Trd&);
  void AddThis         (const G4Trap&);
  void AddThis         (const G4Sphere&);
  void AddThis         (const G4Para&);
  void AddThis         (const G4Torus&);
  void AddThis         (const G4VSolid&);
 
  void BeginPrimitives (const G4Transform3D& objectTransformation);
  void EndPrimitives   ();
  void PreAddThis      (const G4Transform3D& objectTransformation,
			const G4VisAttributes& visAttribs);
  void PostAddThis     ();
  static G4int GetSceneCount ();
  ONode GetRootNode (); 

private:
  void               ClearStore        ();
  void               RequestPrimitives (const G4VSolid& solid);
  G4VGraphicsSystem& fSystem;          // Graphics system for this scene.
  ONode              fRootGoNode;      // Root ONode for this scene.
  static G4int       fSceneIdCount;    // static counter for Wo scenes.
  static G4int       fSceneCount;      // No. of extanct scenes.
  static ONode       fGoNode;          // Current ONode.
  static OColormap   fOColormap;
  void               SetColour         (const G4Colour&);
  char*              nodeName;
};

inline G4int G4GoScene::GetSceneCount () {
  return fSceneCount;
}


#endif

#endif

