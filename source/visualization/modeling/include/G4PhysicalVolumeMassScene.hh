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
// $Id: G4PhysicalVolumeMassScene.hh,v 1.11 2010-05-30 11:23:25 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  12th September 2004

// Class Description:
//
// Calculates the mass of a geometry tree taking into account daughters
// up to the depth specified in the G4PhysicalVolumeModel.  Culling is
// ignored so that all volumes are seen.
//
// Do not use this for a "parallel world" for which materials are not
// defined.  Use only for the material world.
//
// The calculation is quite tricky, since it involves subtracting the
// mass of that part of the mother that is occupied by each daughter and
// then adding the mass of the daughter, and so on down the heirarchy.
//
// Usage for a given G4PhysicalVolumeModel* pvModel:
//   G4PhysicalVolumeMassScene massScene(pvModel);
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
  G4PhysicalVolumeMassScene (G4PhysicalVolumeModel*);
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
  void AddCompound (const G4VDigi&) {}
  void AddCompound (const G4THitsMap<G4double>&) {}

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D&) {}
  virtual void EndPrimitives () {}
  virtual void BeginPrimitives2D (const G4Transform3D&) {}
  virtual void EndPrimitives2D () {}
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
  G4PhysicalVolumeModel* fpPVModel;
  G4double fVolume;
  G4double fMass;
  G4VPhysicalVolume* fpLastPV;
  G4int fPVPCount;
  G4int fLastDepth;
  G4double fLastDensity;
  std::deque<G4double> fDensityStack;
};

#endif
