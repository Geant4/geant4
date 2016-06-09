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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4VoxelLeftBreastROGeometry.hh"

#include "G4SystemOfUnits.hh"
#include "G4HumanDummyLeftBreastSD.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"

G4VoxelLeftBreastROGeometry::G4VoxelLeftBreastROGeometry(G4String aString)
 : G4VReadOutGeometry(aString), ROvoxel_phys(0)
{
}

G4VoxelLeftBreastROGeometry::~G4VoxelLeftBreastROGeometry()
{
}

G4VPhysicalVolume* G4VoxelLeftBreastROGeometry::Build()
{
 
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 
					1., 1.*g/mole, 1.*g/cm3);

  G4double worldSizeX = 2.0*m;
  G4double worldSizeY = 2.0*m;
  G4double worldSizeZ = 2.0*m;

  // world volume of ROGeometry ...
  G4Box *ROWorld = new G4Box("ROWorld",
			     worldSizeX/2.,
			     worldSizeY/2.,
			     worldSizeZ/2.);

  G4LogicalVolume *ROWorldLog = new G4LogicalVolume(ROWorld,
						    dummyMat,
						    "ROWorldLog",
						    0,0,0);

  G4VPhysicalVolume *ROWorldPhys = new G4PVPlacement(0,
						     G4ThreeVector(),
						     "ROWorldPhys",
						     ROWorldLog,
						     0,false,0);
    
  // Breast Mother Volume
  G4double rmin = 0.* cm;
  G4double startPhi = 0.* degree;
  G4double spanningPhi = 180. * degree;
  G4double rmax = 5.5 *cm;
  G4double zz = 5. *cm;
  
  // RO geometry is defined for the inner tube of the breast
 
  G4Tubs* ROinnerLeftBreast = new G4Tubs("innerVoxelLeftBreast",
				   rmin, rmax,
				   zz/2., startPhi, spanningPhi);
  
  G4LogicalVolume* ROinnerLeftBreast_log = new G4LogicalVolume(ROinnerLeftBreast,
							dummyMat,
							"InnerLeftBreast",
							 0, 0, 0);
  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateX(-90.* degree);
  matrix -> rotateY(180.* degree);
  matrix -> rotateZ(-18. * degree);


 G4VPhysicalVolume* ROphysInnerLeftBreast = new G4PVPlacement(matrix,
							  G4ThreeVector(10.*cm, 52.* cm, 8.7 *cm),
      			       "physicalVoxelLeftBreast",
  			       ROinnerLeftBreast_log,
			       ROWorldPhys,
			       false,
			       0);

  // Breast RO Geometry - slices along z axis
  G4double rmin_voxelz = 0.* cm;
  G4double rmax_voxelz = 5.5 * cm;
  G4double zz_voxels = 5. * cm;
  G4double startPhi_voxelz = 0.* degree;
  G4double spanningPhi_voxelz = 180. * degree; 
  G4int nslice = 10;

  G4Tubs* ROvoxelz = new G4Tubs("voxel_z", rmin_voxelz, rmax_voxelz, zz_voxels/(2*nslice),startPhi_voxelz, spanningPhi_voxelz);

  G4LogicalVolume* ROvoxelz_log = new G4LogicalVolume(ROvoxelz, dummyMat, "voxelz_log",0 , 0, 0);
  
  G4VPhysicalVolume* ROvoxelz_phys = new G4PVReplica("LinearArray", ROvoxelz_log, ROphysInnerLeftBreast,
						   kZAxis, nslice, zz_voxels/nslice);

  G4double voxel_height = (zz_voxels/nslice);
  
// Breast RO Geometry - slices along phi angle

  G4int n_voxels = 10;

  G4double voxel_phi = spanningPhi_voxelz/n_voxels;

 
  G4Tubs* ROvoxel = new G4Tubs("voxel", rmin_voxelz, rmax_voxelz, voxel_height/2.,
			     startPhi_voxelz,voxel_phi);
			     
  G4LogicalVolume* ROvoxel_log = new G4LogicalVolume(ROvoxel, dummyMat, "voxel_log", 0,0,0);

  ROvoxel_phys = new G4PVReplica("RZPhiSlices",
 				  ROvoxel_log,
 				 ROvoxelz_phys,
				 kPhi, n_voxels, voxel_phi,- (voxel_phi/2.)); 

  //delete matrix;
						  
  G4HumanDummyLeftBreastSD *dummySD = new G4HumanDummyLeftBreastSD;
  ROvoxel_log -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



