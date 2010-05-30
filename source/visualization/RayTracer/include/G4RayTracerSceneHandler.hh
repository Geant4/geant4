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
// $Id: G4RayTracerSceneHandler.hh,v 1.13 2010-05-30 10:21:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
  void AddCompound(const G4VDigi&){}
  void AddCompound(const G4THitsMap<G4double>&) {}

private:
  static G4int    fSceneIdCount;  // Counter for RayTracer scene handlers.
};

#endif
