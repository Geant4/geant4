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
// History:
//	Pedro Arce  
//
//*******************************************************

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "NestedParamDicomDetectorConstruction.hh"
#include "DicomNestedPhantomParameterisation.hh"

NestedParamDicomDetectorConstruction::NestedParamDicomDetectorConstruction() : DicomDetectorConstruction()
{
}

NestedParamDicomDetectorConstruction::~NestedParamDicomDetectorConstruction()
{
}


//-------------------------------------------------------------
void NestedParamDicomDetectorConstruction::ConstructPhantom()
{
#ifdef G4VERBOSE
  G4cout << "NestedParamDicomDetectorConstruction::ConstructPhantom " << G4endl;
#endif

  //----- Replication of Water Phantom Volume.
  //--- Y Slice
  G4String yRepName("RepY");
  G4VSolid* solYRep =
    new G4Box(yRepName,nVoxelX*voxelHalfDimX,voxelHalfDimY,nVoxelZ*voxelHalfDimZ);
  G4LogicalVolume* logYRep =
    new G4LogicalVolume(solYRep,air,yRepName);
  new G4PVReplica(yRepName,logYRep,container_logic,kYAxis,nVoxelY,voxelHalfDimY*2.);

  //--- X Slice
  G4String xRepName("RepX");
  G4VSolid* solXRep =
    new G4Box(xRepName,voxelHalfDimX,voxelHalfDimY,nVoxelZ*voxelHalfDimZ);
  G4LogicalVolume* logXRep =
    new G4LogicalVolume(solXRep,air,xRepName);
  new G4PVReplica(xRepName,logXRep,logYRep,kXAxis,nVoxelX,voxelHalfDimX*2.);

  //----- Voxel solid and logical volumes
  //--- Z Slice
  G4VSolid* solVoxel = 
    new G4Box("phantom",voxelHalfDimX,voxelHalfDimY,voxelHalfDimZ);
  G4LogicalVolume* logicVoxel = new G4LogicalVolume(solVoxel,air,"phantom");

  //
  // Parameterisation for transformation of voxels.
  //  (voxel size is fixed in this example. 
  //    e.g. nested parameterisation handles material and transfomation of voxels.)
  G4ThreeVector voxelSize(voxelHalfDimX,voxelHalfDimY,voxelHalfDimZ);
  DicomNestedPhantomParameterisation* param
    = new DicomNestedPhantomParameterisation(voxelSize,fMaterials);

  new G4PVParameterised("phantom",    // their name
			logicVoxel, // their logical volume
			logXRep,      // Mother logical volume
			kZAxis,       // Are placed along this axis 
			//			  kUndefined,        // Are placed along this axis 
			nVoxelZ,      // Number of cells
			param);       // Parameterisation.
  
  param->SetMaterialIndices( fMateIDs );
  param->SetNoVoxel( nVoxelX, nVoxelY, nVoxelZ );

}
