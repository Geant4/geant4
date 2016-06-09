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
// $Id: G4RayTracerSceneHandler.hh,v 1.8 2005/06/02 17:43:46 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $

// John Allison  17th March 2000

#ifndef G4RAYTRACERSCENEHANDLER_HH
#define G4RAYTRACERSCENEHANDLER_HH

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
  void AddPrimitive(const G4Scale&){}

  void AddSolid(const G4Box&){}
  void AddSolid(const G4Cons&){}
  void AddSolid(const G4Tubs&){}
  void AddSolid(const G4Trd&){}
  void AddSolid(const G4Trap&){}
  void AddSolid(const G4Sphere&){}
  void AddSolid(const G4Para&){}
  void AddSolid(const G4Torus&){}
  void AddSolid(const G4Polycone&){}
  void AddSolid(const G4Polyhedra&){}
  void AddSolid(const G4VSolid&){}
  void AddCompound(const G4VTrajectory&){}
  void AddCompound(const G4VHit&){}

private:
  static G4int    fSceneIdCount;  // Counter for RayTracer scene handlers.
};

#endif
