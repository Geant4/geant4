// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayTracerSceneHandler.hh,v 1.1 2000-02-23 16:03:52 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// John Allison  17th March 2000

#ifndef RAYTRACERSCENEHANDLER_HH
#define RAYTRACERSCENEHANDLER_HH

#include "G4VSceneHandler.hh"

class G4RayTracerSceneHandler: public G4VSceneHandler {

public:

  G4RayTracerSceneHandler(G4VGraphicsSystem& system,
			  const G4String& name = "");
  virtual ~G4RayTracerSceneHandler();

  void AddPrimitive(const G4Polyline&){}
  void AddPrimitive(const G4Text&){}
  void AddPrimitive(const G4Circle&){}
  void AddPrimitive(const G4Square&){}
  void AddPrimitive(const G4Polyhedron&){}
  void AddPrimitive(const G4NURBS&){}
  void AddPrimitive(const G4Polymarker&){}

  void AddThis(const G4Box&){}
  void AddThis(const G4Cons&){}
  void AddThis(const G4Tubs&){}
  void AddThis(const G4Trd&){}
  void AddThis(const G4Trap&){}
  void AddThis(const G4Sphere&){}
  void AddThis(const G4Para&){}
  void AddThis(const G4Torus&){}
  void AddThis(const G4Polycone&){}
  void AddThis(const G4Polyhedra&){}
  void AddThis(const G4VSolid&){}

  static G4int GetSceneCount();

private:
  static G4int    fSceneIdCount;  // Counter for RayTracer scene handlers.
  static G4int    fSceneCount;    // No. of extanct scenes.
};

#endif
