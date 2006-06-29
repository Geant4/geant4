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
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
//    ************************************
//    *                                  *
//    *      RemSimROGeometry.cc         *
//    *                                  *
//    ************************************
//
// $Id: RemSimROGeometry.cc,v 1.6 2006-06-29 16:24:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "RemSimROGeometry.hh"
#include "RemSimDummySD.hh"
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

RemSimROGeometry::RemSimROGeometry(G4double astronautDimX,
                                   G4double astronautDimY,
                                   G4double astronautDimZ,
                                   G4int numberOfVoxelsZ, 
                                   G4double trans):

  astronautDimensionX(astronautDimX),
  astronautDimensionY(astronautDimY),
  astronautDimensionZ(astronautDimZ),
  numberOfVoxelsAlongZ(numberOfVoxelsZ), 
  translation(trans),
  ROAstronautZDivisionPhys(0)
{
}

RemSimROGeometry::~RemSimROGeometry()
{
}

G4VPhysicalVolume* RemSimROGeometry::Build()
{
  // A dummy material is used to fill the volumes of the readout geometry.
  // (It will be allowed to set a NULL pointer in volumes of such virtual
  // division in future, since this material is irrelevant for tracking.)

  G4Material* dummyMat = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

  G4double worldDimensionX = 30.0*m;
  G4double worldDimensionY = 30.0*m;
  G4double worldDimensionZ = 30.0*m;

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
  
  //  ROGeometry corresponding to the phantom  
  G4Box *ROAstronaut = new G4Box("ROAstronaut", 
			       astronautDimensionX/2.,
			       astronautDimensionY/2., 
			       astronautDimensionZ/2.);

  G4LogicalVolume *ROAstronautLog = new G4LogicalVolume(ROAstronaut,
						      dummyMat,
						      "ROAstronautLog",
						      0,0,0);

  G4VPhysicalVolume *ROAstronautPhys = new G4PVPlacement(0,
						       G4ThreeVector(0.,0.,translation),
						       "AstronautPhys",
                                                       ROAstronautLog,
                                                       ROWorldPhys,
                                                       false,0);
  
  // Z division ...
  G4double voxelSizeZ = astronautDimensionZ/numberOfVoxelsAlongZ;
  
  G4Box *ROAstronautZDivision = new G4Box("ROAstronautZDivision",
					  astronautDimensionX/2,
					  astronautDimensionY/2.,
					  voxelSizeZ/2.);

  G4LogicalVolume *ROAstronautZDivisionLog = new G4LogicalVolume(ROAstronautZDivision,
							       dummyMat,
							       "ROAstronautZDivisionLog",
							       0,0,0);

  ROAstronautZDivisionPhys = new G4PVReplica("ROAstronautXDivisionPhys",
                                               ROAstronautZDivisionLog,
                                               ROAstronautPhys,
                                               kZAxis,
                                               numberOfVoxelsAlongZ,
                                               voxelSizeZ);
  RemSimDummySD *dummySD = new RemSimDummySD;
  ROAstronautZDivisionLog -> SetSensitiveDetector(dummySD);

  return ROWorldPhys;
}



