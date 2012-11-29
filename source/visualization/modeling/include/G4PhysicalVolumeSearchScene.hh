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
// $Id$
//
// 
// John Allison  10th August 1998.
// An artificial scene to find physical volumes.

#ifndef G4PHYSICALVOLUMESEARCHSCENE_HH
#define G4PHYSICALVOLUMESEARCHSCENE_HH

#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"
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
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PhysicalVolumeModel.hh"

class G4PhysicalVolumeSearchScene: public G4VGraphicsScene {

public:
  G4PhysicalVolumeSearchScene
  (G4PhysicalVolumeModel*,
   const G4String& requiredPhysicalVolumeName,
   G4int requiredCopyNo);
  virtual ~G4PhysicalVolumeSearchScene ();
  void PreAddSolid (const G4Transform3D& objectTransformation,
		    const G4VisAttributes&);
  void PostAddSolid ();
  void AddSolid (const G4Box& solid) {FindVolume (solid);}
  void AddSolid (const G4Cons & solid) {FindVolume (solid);}
  void AddSolid (const G4Tubs& solid) {FindVolume (solid);}
  void AddSolid (const G4Trd& solid) {FindVolume (solid);}
  void AddSolid (const G4Trap& solid) {FindVolume (solid);}
  void AddSolid (const G4Sphere& solid) {FindVolume (solid);}
  void AddSolid (const G4Para& solid) {FindVolume (solid);}
  void AddSolid (const G4Torus& solid) {FindVolume (solid);}
  void AddSolid (const G4Polycone& solid) {FindVolume (solid);}
  void AddSolid (const G4Polyhedra& solid) {FindVolume (solid);}
  void AddSolid (const G4VSolid& solid) {FindVolume (solid);}
  void AddCompound (const G4VTrajectory&) {}
  void AddCompound (const G4VHit&) {}
  void AddCompound (const G4VDigi&) {}
  void AddCompound (const G4THitsMap<G4double>&) {}
  
  const std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>&
                        GetFoundFullPVPath     () const;
  G4int                 GetFoundDepth          () const;
  G4VPhysicalVolume*    GetFoundVolume         () const;
  const G4Transform3D&  GetFoundTransformation () const;

  ////////////////////////////////////////////////////////////////
  // Functions not used but required by the abstract interface.

  virtual void BeginPrimitives (const G4Transform3D&) {}
  virtual void EndPrimitives () {}
  virtual void BeginPrimitives2D (const G4Transform3D&) {}
  virtual void EndPrimitives2D () {}
  virtual void AddPrimitive (const G4Polyline&)   {}
  virtual void AddPrimitive (const G4Scale&)       {}
  virtual void AddPrimitive (const G4Text&)       {}
  virtual void AddPrimitive (const G4Circle&)     {}
  virtual void AddPrimitive (const G4Square&)     {}
  virtual void AddPrimitive (const G4Polymarker&) {}
  virtual void AddPrimitive (const G4Polyhedron&) {}
  virtual void AddPrimitive (const G4NURBS&)      {}

private:
  void FindVolume (const G4VSolid&);
  const G4PhysicalVolumeModel* fpPVModel;
  G4String             fRequiredPhysicalVolumeName;
  G4int                fRequiredCopyNo;
  const G4Transform3D* fpCurrentObjectTransformation;
  std::vector<G4PhysicalVolumeModel::G4PhysicalVolumeNodeID>
                       fFoundFullPVPath;            // Found full path.
  G4int                fFoundDepth;                 // Found depth.
  G4VPhysicalVolume*   fpFoundPV;                   // Found physical volume.
  G4LogicalVolume*     fpFoundLV;                   // Found logical volume.
  G4Transform3D        fFoundObjectTransformation;  // Found transformation.
  G4bool               fMultipleOccurrence;
};

#include "G4PhysicalVolumeSearchScene.icc"

#endif
