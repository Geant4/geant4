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
//
// 
// John Allison  15th May 2014
// A base class for "pseudo" scenes/graphics systems.

#ifndef G4PSEUDOSCENE_HH
#define G4PSEUDOSCENE_HH

#include "G4VGraphicsScene.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4TessellatedSolid.hh"

class G4PseudoScene: public G4VGraphicsScene {

public:

  G4PseudoScene(): fpCurrentObjectTransformation(0) {}
  virtual ~G4PseudoScene () {}
  void PreAddSolid (const G4Transform3D& objectTransformation,
		    const G4VisAttributes&)
  {fpCurrentObjectTransformation = &objectTransformation;}
  void PostAddSolid () {}
  // From geometry/solids/CSG
  void AddSolid (const G4Box&    solid) {ProcessVolume (solid);}
  void AddSolid (const G4Cons&   solid) {ProcessVolume (solid);}
  void AddSolid (const G4Orb&    solid) {ProcessVolume (solid);}
  void AddSolid (const G4Para&   solid) {ProcessVolume (solid);}
  void AddSolid (const G4Sphere& solid) {ProcessVolume (solid);}
  void AddSolid (const G4Torus&  solid) {ProcessVolume (solid);}
  void AddSolid (const G4Trap&   solid) {ProcessVolume (solid);}
  void AddSolid (const G4Trd&    solid) {ProcessVolume (solid);}
  void AddSolid (const G4Tubs&   solid) {ProcessVolume (solid);}
  // From geometry/solids/specific
  void AddSolid (const G4Ellipsoid&        solid) {ProcessVolume (solid);}
  void AddSolid (const G4Polycone&         solid) {ProcessVolume (solid);}
  void AddSolid (const G4Polyhedra&        solid) {ProcessVolume (solid);}
  void AddSolid (const G4TessellatedSolid& solid) {ProcessVolume (solid);}
  // For solids not above
  void AddSolid (const G4VSolid& solid) {ProcessVolume (solid);}
  // Compounds
  void AddCompound (const G4VTrajectory&) {}
  void AddCompound (const G4VHit&)        {}
  void AddCompound (const G4VDigi&)       {}
  void AddCompound (const G4THitsMap<G4double>&)     {}
  void AddCompound (const G4THitsMap<G4StatDouble>&) {}
  void AddCompound (const G4Mesh&);  // Catches mesh if special mesh rendering set
  // Primitives
  void BeginPrimitives   (const G4Transform3D&) {}
  void EndPrimitives     ()                     {}
  void BeginPrimitives2D (const G4Transform3D&) {}
  void EndPrimitives2D   ()                     {}
  void AddPrimitive (const G4Polyline&)   {}
  void AddPrimitive (const G4Text&)       {}
  void AddPrimitive (const G4Circle&)     {}
  void AddPrimitive (const G4Square&)     {}
  void AddPrimitive (const G4Polymarker&) {}
  void AddPrimitive (const G4Polyhedron&) {}

  void AddPrimitive (const G4Plotter&)    {}

protected:

  virtual void ProcessVolume (const G4VSolid&);
  const G4Transform3D* fpCurrentObjectTransformation;
};

#endif
