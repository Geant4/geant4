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
// $Id: G4GoSceneHandler.hh,v 1.9 2002/12/11 15:58:00 johna Exp $
// GEANT4 tag $Name: geant4-05-00 $
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
  ////////////////////////////////////////////////////////////////
  // Explicitly invoke base class methods to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive (const G4Scale& scale) {
    G4VSceneHandler::AddPrimitive (scale);
  }

  void AddThis         (const G4Box&);
  void AddThis         (const G4Cons&);
  void AddThis         (const G4Tubs&);
  void AddThis         (const G4Trd&);
  void AddThis         (const G4Trap&);
  void AddThis         (const G4Sphere&);
  void AddThis         (const G4Para&);
  void AddThis         (const G4Torus&);
  void AddThis         (const G4Polycone&);
  void AddThis         (const G4Polyhedra&);
  void AddThis         (const G4VSolid&);
  void AddThis         (const G4VTrajectory&);
  void AddThis         (const G4VHit&);
 
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
  static G4int       fSceneCount;      // No. of extanct scene handlers.
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

