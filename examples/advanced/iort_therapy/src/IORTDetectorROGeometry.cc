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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include "IORTDetectorROGeometry.hh"
#include "IORTDummySD.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

/////////////////////////////////////////////////////////////////////////////
IORTDetectorROGeometry::IORTDetectorROGeometry(G4String aString, G4ThreeVector pos,
								 G4double detectorDimX,
								 G4double detectorDimY,
								 G4double detectorDimZ,
								 G4int numberOfVoxelsX,
								 G4int numberOfVoxelsY,
								 G4int numberOfVoxelsZ):
   
    G4VReadOutGeometry(aString),
    detectorToWorldPosition(pos),
    detectorSizeX(detectorDimX),
    detectorSizeY(detectorDimY),
    detectorSizeZ(detectorDimZ),
    numberOfVoxelsAlongX(numberOfVoxelsX),
    numberOfVoxelsAlongY(numberOfVoxelsY),
    numberOfVoxelsAlongZ(numberOfVoxelsZ)
{
}

/////////////////////////////////////////////////////////////////////////////
IORTDetectorROGeometry::~IORTDetectorROGeometry()
{
}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* IORTDetectorROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  G4double worldSizeX = 400.0 *cm;   // 200.0 *cm;
  G4double worldSizeY = 400.0 *cm;   // 200.0 *cm;
  G4double worldSizeZ = 400.0 *cm;   // 200.0 *cm;

  G4double halfDetectorSizeX = detectorSizeX;
  G4double halfDetectorSizeY = detectorSizeY;
  G4double halfDetectorSizeZ = detectorSizeZ;

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

  // Detector ROGeometry 
  G4Box *RODetector = new G4Box("RODetector", 
			       halfDetectorSizeX, 
			       halfDetectorSizeY, 
			       halfDetectorSizeZ);

  G4LogicalVolume *RODetectorLog = new G4LogicalVolume(RODetector,
						       dummyMat,
						      "RODetectorLog",
						      0,0,0);
  
  G4VPhysicalVolume *RODetectorPhys = new G4PVPlacement(0,
                            detectorToWorldPosition,
							"DetectorPhys",
							RODetectorLog,
							ROWorldPhys,
							false,0);
  
  
  // Division along X axis: the detector is divided in slices along the X axis
  
  G4double halfXVoxelSizeX = halfDetectorSizeX/numberOfVoxelsAlongX;
  G4double halfXVoxelSizeY = halfDetectorSizeY;
  G4double halfXVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelXThickness = 2*halfXVoxelSizeX;

  G4Box *RODetectorXDivision = new G4Box("RODetectorXDivision",
					 halfXVoxelSizeX,
					 halfXVoxelSizeY,
					 halfXVoxelSizeZ);
  
  G4LogicalVolume *RODetectorXDivisionLog = new G4LogicalVolume(RODetectorXDivision,
							       dummyMat,
							       "RODetectorXDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *RODetectorXDivisionPhys = new G4PVReplica("RODetectorXDivisionPhys",
                                                              RODetectorXDivisionLog,
                                                              RODetectorPhys,
                                                              kXAxis,
                                                              numberOfVoxelsAlongX,
                                                              voxelXThickness);

  // Division along Y axis: the slices along the X axis are divided along the Y axis

  G4double halfYVoxelSizeX = halfXVoxelSizeX;
  G4double halfYVoxelSizeY = halfDetectorSizeY/numberOfVoxelsAlongY;
  G4double halfYVoxelSizeZ = halfDetectorSizeZ;
  G4double voxelYThickness = 2*halfYVoxelSizeY;

  G4Box *RODetectorYDivision = new G4Box("RODetectorYDivision",
					halfYVoxelSizeX, 
					halfYVoxelSizeY,
					halfYVoxelSizeZ);

  G4LogicalVolume *RODetectorYDivisionLog = new G4LogicalVolume(RODetectorYDivision,
							       dummyMat,
							       "RODetectorYDivisionLog",
							       0,0,0);
 
  G4VPhysicalVolume *RODetectorYDivisionPhys = new G4PVReplica("RODetectorYDivisionPhys",
							      RODetectorYDivisionLog,
							      RODetectorXDivisionPhys,
							      kYAxis,
							      numberOfVoxelsAlongY,
							      voxelYThickness);
  
  // Division along Z axis: the slices along the Y axis are divided along the Z axis

  G4double halfZVoxelSizeX = halfXVoxelSizeX;
  G4double halfZVoxelSizeY = halfYVoxelSizeY;
  G4double halfZVoxelSizeZ = halfDetectorSizeZ/numberOfVoxelsAlongZ;
  G4double voxelZThickness = 2*halfZVoxelSizeZ;
 
  G4Box *RODetectorZDivision = new G4Box("RODetectorZDivision",
					halfZVoxelSizeX,
					halfZVoxelSizeY, 
					halfZVoxelSizeZ);
 
  G4LogicalVolume *RODetectorZDivisionLog = new G4LogicalVolume(RODetectorZDivision,
							       dummyMat,
							       "RODetectorZDivisionLog",
							       0,0,0);
 
  RODetectorZDivisionPhys = new G4PVReplica("RODetectorZDivisionPhys",
					   RODetectorZDivisionLog,
					   RODetectorYDivisionPhys,
					   kZAxis,
					   numberOfVoxelsAlongZ,
					   voxelZThickness);

  IORTDummySD *dummySD = new IORTDummySD;
  RODetectorZDivisionLog -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



