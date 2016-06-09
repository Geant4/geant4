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
// $Id: G4PhysicalVolumeMassScene.hh,v 1.4 2005/01/27 20:06:35 johna Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// John Allison  12th September 2004

// Class Description:
//
// Calculates the mass of a geometry tree taking into account daughters
// up to the depth specified in the G4PhysicalVolumeModel.  Culling is
// ignored so that all volumes are seen.
//
// The calculation is quite tricky, since it involves subtracting the
// mass of that part of the mother that is occupied by each daughter and
// then adding the mass of the daughter, and so on down the heirarchy.
//
// Usage for a given G4PhysicalVolumeModel* pvModel:
//   G4PhysicalVolumeMassScene massScene;
//   massScene.EstablishSpecials (*pvModel);
//   pvModel->DescribeYourselfTo (massScene);
//   G4double volume = massScene.GetVolume();
//   G4double mass = massScene.GetMass();
//   massScene.Reset();
// See, for example, G4ASCIITreeSceneHandler::EndModeling().

#ifndef G4PHYSICALVOLUMEMASSSCENE_HH
#define G4PHYSICALVOLUMEMASSSCENE_HH

#include "G4VGraphicsScene.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include <deque>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4PhysicalVolumeModel;
class G4Material;

class G4PhysicalVolumeMassScene: public G4VGraphicsScene {

public:
  G4PhysicalVolumeMassScene ();
  virtual ~G4PhysicalVolumeMassScene ();

public: // With description

  G4double GetVolume () const {return fVolume;}
  // Overall volume.

  G4double GetMass () const {return fMass;}
  // Mass of whole tree, i.e., accounting for all daughters.

  void Reset ();
  // Reset for subsequent re-use.

public:

  // Force execution of AccrueMass for all solids...
  void PreAddSolid (const G4Transform3D&, const G4VisAttributes&) {}
  void PostAddSolid () {}
  void AddSolid (const G4Box& s) {AccrueMass (s);}
  void AddSolid (const G4Cons & s) {AccrueMass (s);}
  void AddSolid (const G4Tubs& s) {AccrueMass (s);}
  void AddSolid (const G4Trd& s) {AccrueMass (s);}
  void AddSolid (const G4Trap& s) {AccrueMass (s);}
  void AddSolid (const G4Sphere& s) {AccrueMass (s);}
  void AddSolid (const G4Para& s) {AccrueMass (s);}
  void AddSolid (const G4Torus& s) {AccrueMass (s);}
  void AddSolid (const G4Polycone& s) {AccrueMass (s);}
  void AddSolid (const G4Polyhedra& s) {AccrueMass (s);}
  void AddSolid (const G4VSolid& s) {AccrueMass (s);}
  void AddCompound (const G4VTrajectory&) {}
  void AddCompound (const G4VHit&) {}
  void EstablishSpecials (G4PhysicalVolumeModel&);

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D&) {}
  virtual void EndPrimitives () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Scale&)      {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}

private:
  void AccrueMass (const G4VSolid&);
  G4double fVolume;
  G4double fMass;
  G4VPhysicalVolume* fpLastPV;
  G4int fPVPCount;
  G4int fLastDepth;
  G4double fLastDensity;
  std::deque<G4double> fDensityStack;
  G4int                fCurrentDepth;  // Current depth of geom. hierarchy.
  G4VPhysicalVolume*   fpCurrentPV;    // Current physical volume.
  G4LogicalVolume*     fpCurrentLV;    // Current logical volume.
  G4Material*      fpCurrentMaterial;  // Current material.
};

#endif
