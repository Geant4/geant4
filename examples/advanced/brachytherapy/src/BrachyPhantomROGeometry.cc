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
// $Id: BrachyPhantomROGeometry.cc,v 1.10 2006-05-12 13:23:48 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyPhantomROGeometry.hh"
#include "BrachyDummySD.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

BrachyPhantomROGeometry::BrachyPhantomROGeometry(G4String aString,
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

BrachyPhantomROGeometry::~BrachyPhantomROGeometry()
{
}

G4VPhysicalVolume* BrachyPhantomROGeometry::Build()
{
 
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 
					1., 1.*g/mole, 1.*g/cm3);

  G4double worldSizeX = 4.0*m;
  G4double worldSizeY = 4.0*m;
  G4double worldSizeZ = 4.0*m;

  G4double halfPhantomSizeX = phantomSizeX;
  G4double halfPhantomSizeZ = phantomSizeZ;
  G4double halfPhantomSizeY = phantomSizeY;

  // variables for x division ...
  G4double halfVoxelSize = halfPhantomSizeX/numberOfVoxelsAlongX;
  G4double halfSize = halfPhantomSizeZ;
  G4double voxelThickness = 2*halfVoxelSize;

  // world volume of ROGeometry ...
  G4Box *ROWorld = new G4Box("ROWorld",
			     worldSizeX,
			     worldSizeY,
			     worldSizeZ);

  G4LogicalVolume *ROWorldLog = new G4LogicalVolume(ROWorld,
						    dummyMat,
						    "ROWorldLog",
						    0,0,0);

  G4VPhysicalVolume *ROWorldPhys = new G4PVPlacement(0,
						     G4ThreeVector(),
						     "ROWorldPhys",
						     ROWorldLog,
						     0,false,0);
  
  // Phantom ROGeometry ... 
  G4Box *ROPhantom = new G4Box("ROPhantom", 
			       halfPhantomSizeX, 
			       halfPhantomSizeY, 
			       halfPhantomSizeZ);

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
 
  // ROGeometry: Voxel division
 
  // X division first... 
  G4Box *ROPhantomXDivision = new G4Box("ROPhantomXDivision",
					halfVoxelSize,
					halfSize,
					halfSize);

  G4LogicalVolume *ROPhantomXDivisionLog = new G4LogicalVolume(ROPhantomXDivision,
							       dummyMat,
							       "ROPhantomXDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomXDivisionPhys = new G4PVReplica("ROPhantomXDivisionPhys",
                                                              ROPhantomXDivisionLog,
                                                              ROPhantomPhys,
                                                              kXAxis,
                                                              numberOfVoxelsAlongX,
                                                              voxelThickness);
  // ...then Z division
  
  G4Box *ROPhantomZDivision = new G4Box("ROPhantomZDivision",
					halfVoxelSize,
					halfSize, 
					halfVoxelSize);

  G4LogicalVolume *ROPhantomZDivisionLog = new G4LogicalVolume(ROPhantomZDivision,
							       dummyMat,
							       "ROPhantomZDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomZDivisionPhys = new G4PVReplica("ROPhantomZDivisionPhys",
							      ROPhantomZDivisionLog,
							      ROPhantomXDivisionPhys,
							      kZAxis,
							      numberOfVoxelsAlongZ,
							      voxelThickness);
  // ...then Y  division

  G4Box *ROPhantomYDivision = new G4Box("ROPhantomYDivision",
					halfVoxelSize, 
					halfVoxelSize,
					halfVoxelSize);

  G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
							       dummyMat,
							       "ROPhantomYDivisionLog",
							       0,0,0);
 
  ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
							      ROPhantomYDivisionLog,
							      ROPhantomZDivisionPhys,
							      kYAxis,
							      numberOfVoxelsAlongY,
							      voxelThickness);
  BrachyDummySD *dummySD = new BrachyDummySD;
  ROPhantomYDivisionLog -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



