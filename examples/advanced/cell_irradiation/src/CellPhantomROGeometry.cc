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
//
//    ************************************
//    *                                  *
//    *     CellPhantomROGeometry.cc     *
//    *                                  *
//    ************************************
//
//
//
#include "CellPhantomROGeometry.hh"
#include "CellDummySD.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

CellPhantomROGeometry::CellPhantomROGeometry(G4String aString,
                                               G4double phantomDimX,
					       G4double phantomDimY,
                                               G4double phantomDimZ,
                                               G4int numberOfVoxelsX,
                                               G4int numberOfVoxelsY,
					       G4int numberOfVoxelsZ):
  G4VReadOutGeometry(aString),
  phantomDimensionX(phantomDimX),
  phantomDimensionY(phantomDimY),
  phantomDimensionZ(phantomDimZ),
  numberOfVoxelsAlongX(numberOfVoxelsX),
  numberOfVoxelsAlongY(numberOfVoxelsY),
  numberOfVoxelsAlongZ(numberOfVoxelsZ)
{
}

CellPhantomROGeometry::~CellPhantomROGeometry()
{
}

G4VPhysicalVolume* CellPhantomROGeometry::Build()
{
 
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  G4double worldDimensionX = 20.0*m;
  G4double worldDimensionY = 20.0*m;
  G4double worldDimensionZ = 20.0*m;

  G4double halfPhantomDimensionX = phantomDimensionX/2.;
  G4double halfPhantomDimensionZ = phantomDimensionZ/2.;
  G4double halfPhantomDimensionY = phantomDimensionX/2.;

  // variables for x division ...
  G4double halfXVoxelDimensionX = halfPhantomDimensionX/numberOfVoxelsAlongX;
  G4double halfXVoxelDimensionZ = halfPhantomDimensionZ;
  G4double halfXVoxelDimensionY = halfPhantomDimensionY;
  G4double voxelXThickness = 2*halfXVoxelDimensionX;

 // variables for z division ...
  G4double halfZVoxelDimensionX = halfPhantomDimensionX/numberOfVoxelsAlongX;
  G4double halfZVoxelDimensionZ = halfPhantomDimensionZ/numberOfVoxelsAlongZ;
  G4double halfZVoxelDimensionY = halfPhantomDimensionY;
  G4double voxelZThickness = 2*halfZVoxelDimensionZ;


 // variables for y division ...
  G4double halfYVoxelDimensionX = halfPhantomDimensionX/numberOfVoxelsAlongX;
  G4double halfYVoxelDimensionZ = halfPhantomDimensionZ/numberOfVoxelsAlongZ;
  G4double halfYVoxelDimensionY = halfPhantomDimensionY/numberOfVoxelsAlongY;
  G4double voxelYThickness = 2*halfYVoxelDimensionY;


  // world volume of ROGeometry ...
  G4Box *ROWorld = new G4Box("ROWorld",
			     worldDimensionX/2,
			     worldDimensionY/2,
			     worldDimensionZ/2);

  G4LogicalVolume *ROWorldLog = new G4LogicalVolume(ROWorld,
						    dummyMat,
						    "ROWorldLog",
						    0,0,0);

  G4VPhysicalVolume *ROWorldPhys = new G4PVPlacement(0,
						     G4ThreeVector(),
						     "ROWorldPhys",
						     ROWorldLog,
						     0,false,0);
  
  // phantom ROGeometry ... 
  G4Box *ROPhantom = new G4Box("ROPhantom", 
			       halfPhantomDimensionX, 
			       halfPhantomDimensionY, 
			       halfPhantomDimensionZ);

  G4LogicalVolume *ROPhantomLog = new G4LogicalVolume(ROPhantom,
						      dummyMat,
						      "ROPhantomLog",
						      0,0,0);

  G4VPhysicalVolume *ROPhantomPhys = new G4PVPlacement(0,
						       G4ThreeVector(),
						       "PhantomPhys",
                                                       ROPhantomLog,
                                                       ROWorldPhys,
                                                       false,0);
 
  // ROGeomtry: Voxel division
 
  // X division first... 
  G4Box *ROPhantomXDivision = new G4Box("ROPhantomXDivision",
					halfXVoxelDimensionX,
					halfXVoxelDimensionY,
					halfXVoxelDimensionZ);

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
  // ...then Z division
  
  G4Box *ROPhantomZDivision = new G4Box("ROPhantomZDivision",
					halfZVoxelDimensionX,
					halfZVoxelDimensionY, 
					halfZVoxelDimensionZ);

  G4LogicalVolume *ROPhantomZDivisionLog = new G4LogicalVolume(ROPhantomZDivision,
							       dummyMat,
							       "ROPhantomZDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomZDivisionPhys = new G4PVReplica("ROPhantomZDivisionPhys",
							      ROPhantomZDivisionLog,
							      ROPhantomXDivisionPhys,
							      kZAxis,
							      numberOfVoxelsAlongZ,
							      voxelZThickness);
  // ...then Y  division

  G4Box *ROPhantomYDivision = new G4Box("ROPhantomYDivision",
					halfYVoxelDimensionX, 
					halfYVoxelDimensionY,
					halfYVoxelDimensionZ);

  G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
							       dummyMat,
							       "ROPhantomYDivisionLog",
							       0,0,0);
 
  ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
							      ROPhantomYDivisionLog,
							      ROPhantomZDivisionPhys,
							      kYAxis,
							      numberOfVoxelsAlongY,
							      voxelYThickness);
  CellDummySD *dummySD = new CellDummySD;
  ROPhantomYDivisionLog -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



