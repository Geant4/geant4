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
// $Id: HadrontherapyPhantomROGeometry.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyPhantomROGeometry.hh"
#include "HadrontherapyDummySD.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

HadrontherapyPhantomROGeometry::HadrontherapyPhantomROGeometry(G4String aString,
							       G4double phantomDimX,
							       G4double phantomDimY,
							       G4double phantomDimZ,
							       G4int numberOfVoxelsX,
							       G4int numberOfVoxelsY,
							       G4int numberOfVoxelsZ):
  G4VReadOutGeometry(aString),
  phantomSizeX(phantomDimX),
  phantomSizeY(phantomDimY),
  phantomSizeZ(phantomDimZ),
  numberOfVoxelsAlongX(numberOfVoxelsX),
  numberOfVoxelsAlongY(numberOfVoxelsY),
  numberOfVoxelsAlongZ(numberOfVoxelsZ)
{
}

HadrontherapyPhantomROGeometry::~HadrontherapyPhantomROGeometry()
{
}

G4VPhysicalVolume* HadrontherapyPhantomROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  G4double worldSizeX = 200.0 *cm;
  G4double worldSizeY = 200.0 *cm;
  G4double worldSizeZ = 200.0 *cm;

  G4double halfPhantomSizeX = phantomSizeX;
  G4double halfPhantomSizeY = phantomSizeY;
  G4double halfPhantomSizeZ = phantomSizeZ;

  // World volume of ROGeometry ...
  G4Box* ROWorld = new G4Box("ROWorld",
			     worldSizeX,
			     worldSizeY,
			     worldSizeZ);

  G4LogicalVolume* ROWorldLog = new G4LogicalVolume(ROWorld, dummyMat, 
						    "ROWorldLog", 0,0,0);

  G4VPhysicalVolume* ROWorldPhys = new G4PVPlacement(0,G4ThreeVector(), 
						     "ROWorldPhys", 
						     ROWorldLog, 
						     0,false,0);

  // Phantom ROGeometry 
  G4Box *ROPhantom = new G4Box("ROPhantom", 
			       halfPhantomSizeX, 
			       halfPhantomSizeY, 
			       halfPhantomSizeZ);

  G4LogicalVolume *ROPhantomLog = new G4LogicalVolume(ROPhantom,
						      dummyMat,
						      "ROPhantomLog",
						      0,0,0);
  
  G4VPhysicalVolume *ROPhantomPhys = new G4PVPlacement(0,
						       G4ThreeVector(-180.0 *mm,
								     0.0 *mm, 
								     0.0 *mm),
						       "PhantomPhys",
                                                       ROPhantomLog,
                                                       ROWorldPhys,
                                                       false,0);


  // ROGeomtry: the phantom is divided in voxels along the axis X, Y, Z

  // Division along X axis: the phantom is devided in slices along the X axis

  G4double halfXVoxelSizeX = halfPhantomSizeX/numberOfVoxelsAlongX;
  G4double halfXVoxelSizeY = halfPhantomSizeY;
  G4double halfXVoxelSizeZ = halfPhantomSizeZ;
  G4double voxelXThickness = 2*halfXVoxelSizeX;

  G4Box *ROPhantomXDivision = new G4Box("ROPhantomXDivision",
					halfXVoxelSizeX,
					halfXVoxelSizeY,
					halfXVoxelSizeZ);

  G4LogicalVolume *ROPhantomXDivisionLog = new G4LogicalVolume(ROPhantomXDivision,
							       dummyMat,
							       "ROPhantomXDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomXDivisionPhys = new G4PVReplica("ROPhantomXDivisionPhys",
                                                              ROPhantomXDivisionLog,
                                                              ROPhantomPhys,
                                                              kXAxis,
                                                              numberOfVoxelsAlongX,
                                                              voxelXThickness);

  // Division along Y axis: the slices along the X axis are devided along the Y axis

  G4double halfYVoxelSizeX =  halfXVoxelSizeX;
  G4double halfYVoxelSizeY = halfPhantomSizeY/numberOfVoxelsAlongY;
  G4double halfYVoxelSizeZ = halfPhantomSizeZ;
  G4double voxelYThickness = 2*halfYVoxelSizeY;

  G4Box *ROPhantomYDivision = new G4Box("ROPhantomYDivision",
					halfYVoxelSizeX, 
					halfYVoxelSizeY,
					halfYVoxelSizeZ);

  G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
							       dummyMat,
							       "ROPhantomYDivisionLog",
							       0,0,0);
 
  G4VPhysicalVolume *ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
							      ROPhantomYDivisionLog,
							      ROPhantomXDivisionPhys,
							      kYAxis,
							      numberOfVoxelsAlongY,
							      voxelYThickness);
  
  // Division along Z axis: the slices along the Y axis are devided along the Z axis

  G4double halfZVoxelSizeX = halfXVoxelSizeX;
  G4double halfZVoxelSizeY = halfYVoxelSizeY;
  G4double halfZVoxelSizeZ = halfPhantomSizeZ/numberOfVoxelsAlongZ;
  G4double voxelZThickness = 2*halfZVoxelSizeZ;
 
  G4Box *ROPhantomZDivision = new G4Box("ROPhantomZDivision",
					halfZVoxelSizeX,
					halfZVoxelSizeY, 
					halfZVoxelSizeZ);
 
  G4LogicalVolume *ROPhantomZDivisionLog = new G4LogicalVolume(ROPhantomZDivision,
							       dummyMat,
							       "ROPhantomZDivisionLog",
							       0,0,0);
 
  ROPhantomZDivisionPhys = new G4PVReplica("ROPhantomZDivisionPhys",
					   ROPhantomZDivisionLog,
					   ROPhantomYDivisionPhys,
					   kZAxis,
					   numberOfVoxelsAlongZ,
					   voxelZThickness);

  HadrontherapyDummySD *dummySD = new HadrontherapyDummySD;
  ROPhantomZDivisionLog -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



