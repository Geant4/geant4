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
// $Id: G4LogicalVolumeModel.cc 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// John Allison  26th July 1999.
// Model for logical volumes.

#include "G4LogicalVolumeModel.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
#include "G4DrawVoxels.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VReadOutGeometry.hh"

G4LogicalVolumeModel::G4LogicalVolumeModel
(G4LogicalVolume*            pLV,
 G4int                       soughtDepth,
 G4bool                      booleans,
 G4bool                      voxels,
 G4bool                      readout,
 const G4Transform3D&        modelTransformation,
 const G4ModelingParameters* pMP):
  // Instantiate a G4PhysicalVolumeModel with a G4PVPlacement to
  // represent this logical volume.  It has no rotation and a null
  // translation so that the logical volume will be seen in its own
  // reference system.  It will be added to the physical volume store
  // but it will not be part of the normal geometry heirarchy so it
  // has no mother.
  G4PhysicalVolumeModel
(new G4PVPlacement (0,                   // No rotation.
		    G4ThreeVector(),     // Null traslation.
		    "PhysVol representation of LogVol " + pLV -> GetName (),
		    pLV,
		    0,                   // No mother.
		    false,               // Not "MANY".
		    0),                  // Copy number.
 soughtDepth,
 modelTransformation,
 pMP,
 true),                                  // Use full extent.
  fpLV (pLV),
  fBooleans (booleans),
  fVoxels (voxels),
  fReadout (readout)
{
  fType = "G4LogicalVolumeModel";
  fGlobalTag = fpLV -> GetName ();
  fGlobalDescription = "G4LogicalVolumeModel " + fGlobalTag;
}

G4LogicalVolumeModel::~G4LogicalVolumeModel () {}

void G4LogicalVolumeModel::DescribeYourselfTo
(G4VGraphicsScene& sceneHandler) {

  // Store current modeling parameters and ensure nothing is culled.
  const G4ModelingParameters* tmpMP = fpMP;
  G4ModelingParameters nonCulledMP;
  if (fpMP) nonCulledMP = *fpMP;
  nonCulledMP.SetCulling (false);
  fpMP = &nonCulledMP;    
  G4PhysicalVolumeModel::DescribeYourselfTo (sceneHandler);
  fpMP = tmpMP;

  if (fVoxels) {
    if (fpTopPV->GetLogicalVolume()->GetVoxelHeader()) {
      // Add Voxels.
      G4DrawVoxels dv;
      G4PlacedPolyhedronList* pPPL =
	dv.CreatePlacedPolyhedra (fpTopPV -> GetLogicalVolume ());
      for (size_t i = 0; i < pPPL -> size (); i++) {
	const G4Transform3D& transform = (*pPPL)[i].GetTransform ();
	const G4Polyhedron& polyhedron = (*pPPL)[i].GetPolyhedron ();
	sceneHandler.BeginPrimitives (transform);
	sceneHandler.AddPrimitive (polyhedron);
	sceneHandler.EndPrimitives ();
      }
      delete pPPL;
    }
  }

  if (fReadout) {
    // Draw readout geometry...
    G4VSensitiveDetector* sd = fpLV->GetSensitiveDetector();
    if (sd) {
      G4VReadOutGeometry* roGeom = sd->GetROgeometry();
      if (roGeom) {
	G4VPhysicalVolume* roWorld = roGeom->GetROWorld();
	G4cout << "Readout geometry \"" << roGeom->GetName()
	       << "\" with top physical volume \""
	       << roWorld->GetName()
	       << "\"" << G4endl;
	G4PhysicalVolumeModel pvModel(roWorld);
	pvModel.SetModelingParameters(fpMP);
	pvModel.DescribeYourselfTo(sceneHandler);
      }
    }
  }
}

// This called from G4PhysicalVolumeModel::DescribeAndDescend by the
// virtual function mechanism.
void G4LogicalVolumeModel::DescribeSolid
(const G4Transform3D& theAT,
 G4VSolid* pSol,
 const G4VisAttributes* pVisAttribs,
 G4VGraphicsScene& sceneHandler) {

  if (fBooleans) {
    // Look for "constituents".  Could be a Boolean solid.
    G4VSolid* pSol0 = pSol -> GetConstituentSolid (0);
    if (pSol0) {  // Composite solid...
      G4VSolid* pSol1 = pSol -> GetConstituentSolid (1);
      if (!pSol1) {
	G4Exception
	  ("G4PhysicalVolumeModel::DescribeSolid",
	   "modeling0001", FatalException,
	   "2nd component solid in Boolean is missing.");
      }
      // Draw these constituents white and "forced wireframe"...
      G4VisAttributes constituentAttributes;
      constituentAttributes.SetForceWireframe(true);
      DescribeSolid (theAT, pSol0, &constituentAttributes, sceneHandler);
      DescribeSolid (theAT, pSol1, &constituentAttributes, sceneHandler);
    }
  }

  // In any case draw the original/resultant solid...
  sceneHandler.PreAddSolid (theAT, *pVisAttribs);
  pSol -> DescribeYourselfTo (sceneHandler);
  sceneHandler.PostAddSolid ();
}
