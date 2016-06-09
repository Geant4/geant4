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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ************************************
//    *                                  *
//    *   BrachyPhantomROGeometry.cc    *
//    *                                  *
//    ************************************
//
// $Id: BrachyPhantomROGeometry.cc,v 1.7 2003/11/18 13:31:46 guatelli Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
#include "BrachyPhantomROGeometry.hh"
#include "BrachyDummySD.hh"

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

BrachyPhantomROGeometry::BrachyPhantomROGeometry(G4String aString,
                                                 G4double phantomDimX,
                                                 G4double phantomDimZ,
                                                 G4int numberOfVoxelsX,
                                                 G4int numberOfVoxelsZ):
  G4VReadOutGeometry(aString),
  phantomDimensionX(phantomDimX),
  phantomDimensionZ(phantomDimZ),
  numberOfVoxelsAlongX(numberOfVoxelsX),
  numberOfVoxelsAlongZ(numberOfVoxelsZ)
{
}

BrachyPhantomROGeometry::~BrachyPhantomROGeometry()
{
}

G4VPhysicalVolume* BrachyPhantomROGeometry::Build()
{
 
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  G4double worldDimensionX = 4.0*m;
  G4double worldDimensionY = 4.0*m;
  G4double worldDimensionZ = 4.0*m;

  G4double halfPhantomDimensionX = phantomDimensionX/2;
  G4double halfPhantomDimensionZ = phantomDimensionZ/2;
  G4double halfPhantomDimensionY = phantomDimensionX/2 ;

  // variables for x division ...
  G4double halfXVoxelDimensionX = halfPhantomDimensionX/numberOfVoxelsAlongX;
  G4double halfXVoxelDimensionZ = halfPhantomDimensionZ;
  G4double voxelXThickness = 2*halfXVoxelDimensionX;

  // variables for z and y division ...
  G4double halfZVoxelDimensionX = halfPhantomDimensionX; 
  G4double halfZVoxelDimensionZ = halfPhantomDimensionZ/numberOfVoxelsAlongZ;
  G4double voxelZThickness = 2*halfZVoxelDimensionZ;

  // world volume of ROGeometry ...
  G4Box *ROWorld = new G4Box("ROWorld",
			     worldDimensionX,
			     worldDimensionY,
			     worldDimensionZ);

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
					halfXVoxelDimensionZ,
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
					halfZVoxelDimensionX, 
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
					halfZVoxelDimensionZ, 
					halfZVoxelDimensionZ,
					halfZVoxelDimensionZ);

  G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
							       dummyMat,
							       "ROPhantomYDivisionLog",
							       0,0,0);
 
  ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
							      ROPhantomYDivisionLog,
							      ROPhantomZDivisionPhys,
							      kYAxis,
							      numberOfVoxelsAlongZ,
							      voxelZThickness);
  BrachyDummySD *dummySD = new BrachyDummySD;
  ROPhantomYDivisionLog->SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



