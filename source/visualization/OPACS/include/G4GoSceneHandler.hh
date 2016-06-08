// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GoSceneHandler.hh,v 1.3 1999/05/25 12:04:43 barrand Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
// 
// Guy Barrand 04 November 1996
// Wo scene handler - creates Wo Display lists.

#ifndef G4GOSCENEHANDLER_HH
#define G4GOSCENEHANDLER_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

//Go
#include <ONode.h>
#include <OColormap.h>
//G4
#include "globals.hh"
#include "G4VGraphicsSystem.hh"
#include "G4VSceneHandler.hh"

class G4GoSceneHandler: public G4VSceneHandler {

public:
  G4GoSceneHandler            (G4VGraphicsSystem& system,
			const G4String& name = "");
  virtual ~G4GoSceneHandler   ();
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
  void               ClearTransientStore();
  void               RequestPrimitives (const G4VSolid& solid);
  G4VGraphicsSystem& fSystem;          // Graphics system for this scene.
  ONode              fRootGoNode;            // Root ONode for this scene.
  ONode              fStaticRootGoNode;      // For detector.
  ONode              fTransientRootGoNode;   // For event.
  static G4int       fSceneIdCount;    // static counter for Wo scenes.
  static G4int       fSceneCount;      // No. of extanct scenes.
  static ONode       fGoNode;          // Current ONode.
  static OColormap   fOColormap;
  void               SetColour         (const G4Colour&);
  char*              nodeName;
};

inline G4int G4GoSceneHandler::GetSceneCount () {
  return fSceneCount;
}


#endif

#endif

