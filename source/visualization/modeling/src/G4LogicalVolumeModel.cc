// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LogicalVolumeModel.cc,v 1.2 1999-12-15 14:54:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th July 1999.
// Model for logical volumes.

#include "G4LogicalVolumeModel.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ModelingParameters.hh"
#include "G4VGraphicsScene.hh"
////////////////#include "G4DrawVoxels.hh"

G4LogicalVolumeModel::G4LogicalVolumeModel
(G4LogicalVolume*            pLV,
 G4int                       soughtDepth,
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
		    "PhysVol representaion of LogVol " + pLV -> GetName (),
		    pLV,
		    0,                   // No mother.
		    false,               // Not "MANY".
		    0),                  // Copy number.
 soughtDepth,
 modelTransformation,
 pMP,
 true),                                  // Use full extent.
  fpLV (pLV)
{
  fGlobalTag = fpLV -> GetName ();
  fGlobalDescription = "G4LogicalVolumeModel " + fGlobalTag;
}

G4LogicalVolumeModel::~G4LogicalVolumeModel () {}

G4String G4LogicalVolumeModel::GetCurrentDescription () const {
  return "G4LogicalVolumeModel " + GetCurrentTag ();
}

void G4LogicalVolumeModel::DescribeYourselfTo
(G4VGraphicsScene& sceneHandler) {

  // Store current modeling parameters and ensure nothing is culled.
  const G4ModelingParameters* tpMP = fpMP;
  G4ModelingParameters nonCulledMP (*fpMP);
  nonCulledMP.SetCulling (false);
  fpMP = &nonCulledMP;    

  G4PhysicalVolumeModel::DescribeYourselfTo (sceneHandler);

  fpMP = tpMP;

  /*****************************
  // Add Voxels.
  G4DrawVoxels dv;
  G4PlacedPolyhedronList* pPPL =
    dv.CreatePlacedPolyhedra (fpTopPV -> GetLogicalVolume ());
  for (int i = 0; i < pPPL -> entries (); i++) {
    const G4Transform3D& transform = (*pPPL)[i].GetTransform ();
    const G4Polyhedron& polyhedron = (*pPPL)[i].GetPolyhedron ();
    sceneHandler.BeginPrimitives (transform);
    sceneHandler.AddPrimitive (polyhedron);
    sceneHandler.EndPrimitives ();
  }
  delete pPPL;
  *********************************/
}
