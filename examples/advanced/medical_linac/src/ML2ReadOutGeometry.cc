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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//

#include "ML2ReadOutGeometry.hh"

#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "ML2DummySD.hh"
#include "G4PVReplica.hh"

CML2ReadOutGeometry::CML2ReadOutGeometry() : ROPhyVol(0)
{
    // Build the world volume
    G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4ThreeVector halfSizeWorld, ctr;
    ctr.set(0.*mm, 0.*mm, 0.*mm);
    halfSizeWorld.set(3000.*mm, 3000*mm, 3000*mm);
    G4Box *ROphmWorldB = new G4Box("ROphmWorldG", halfSizeWorld.getX(), halfSizeWorld.getY(), halfSizeWorld.getZ());
    G4LogicalVolume *ROphmWorldLV = new G4LogicalVolume(ROphmWorldB, Vacuum, "ROphmWorldL", 0, 0, 0);
    ROPhyVol= new G4PVPlacement(0, ctr, "ROphmWorldPV", ROphmWorldLV, 0, false, 0);
}

CML2ReadOutGeometry::~CML2ReadOutGeometry(void)
{
    delete ROPhyVol;
}

void CML2ReadOutGeometry::setBuildData(G4ThreeVector ctr, G4ThreeVector hSiz, G4int NVX, G4int NVY, G4int NVZ)
{
    centre=ctr;
    halfSize=hSiz;
    NumberOfVoxelsAlongX=NVX;
    NumberOfVoxelsAlongY=NVY;
    NumberOfVoxelsAlongZ=NVZ;
}

G4VPhysicalVolume* CML2ReadOutGeometry::Build()
{
    // Build RO Zone
    G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Box *ROBox = new G4Box("ROBox", halfSize.getX(), halfSize.getY(), halfSize.getZ());
    G4LogicalVolume *ROLV = new G4LogicalVolume(ROBox, Vacuum, "ROLV", 0, 0, 0);
    G4VPhysicalVolume *ROPV;
    ROPV = new G4PVPlacement(0, centre, "ROPV", ROLV, ROPhyVol, false, 0);

    // ROGeomtry: Voxel division

    G4double halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ;

    halfXVoxelDimensionX=halfSize.getX()/NumberOfVoxelsAlongX;
    halfXVoxelDimensionY=halfSize.getY()/NumberOfVoxelsAlongY;
    halfXVoxelDimensionZ=halfSize.getZ()/NumberOfVoxelsAlongZ;

    G4double voxelXThicknessX = 2*halfXVoxelDimensionX;
    G4double voxelXThicknessY = 2*halfXVoxelDimensionY;
    G4double voxelXThicknessZ = 2*halfXVoxelDimensionZ;


    // X division first... slice along X axis
    G4Box *ROPhantomXDivision = new G4Box("ROPhantomXDivision",
                                          halfXVoxelDimensionX,
                                          halfSize.getY(),
                                          halfSize.getZ());

    G4LogicalVolume *ROPhantomXDivisionLog = new G4LogicalVolume(ROPhantomXDivision,
                                                                 Vacuum,
                                                                 "ROPhantomXDivisionLog",
                                                                 0,0,0);
    G4VPhysicalVolume *ROPhantomXDivisionPhys;
    ROPhantomXDivisionPhys = new G4PVReplica("ROPhantomXDivisionPhys",
                                             ROPhantomXDivisionLog,
                                             ROPV,
                                             kXAxis,
                                             NumberOfVoxelsAlongX,
                                             voxelXThicknessX,
                                             -halfSize.getX());
    // ...then Z division

    G4Box *ROPhantomZDivision = new G4Box("ROPhantomZDivision",
                                          halfXVoxelDimensionX,
                                          halfSize.getY(),
                                          halfXVoxelDimensionZ);

    G4LogicalVolume *ROPhantomZDivisionLog = new G4LogicalVolume(ROPhantomZDivision,
                                                                 Vacuum,
                                                                 "ROPhantomZDivisionLog",
                                                                 0,0,0);
    G4VPhysicalVolume *ROPhantomZDivisionPhys;
    ROPhantomZDivisionPhys = new G4PVReplica("ROPhantomZDivisionPhys",
                                             ROPhantomZDivisionLog,
                                             ROPhantomXDivisionPhys,
                                             kZAxis,
                                             NumberOfVoxelsAlongZ,
                                             voxelXThicknessZ,
                                             -halfSize.getZ());
    // ...then Y  division

    G4Box *ROPhantomYDivision = new G4Box("ROPhantomYDivision",
                                          halfXVoxelDimensionX,
                                          halfXVoxelDimensionY,
                                          halfXVoxelDimensionZ);

    G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
                                                                 Vacuum,
                                                                 "ROPhantomYDivisionLog",
                                                                 0,0,0);
    ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
                                             ROPhantomYDivisionLog,
                                             ROPhantomZDivisionPhys,
                                             kYAxis,
                                             NumberOfVoxelsAlongY,
                                             voxelXThicknessY,
                                             -halfSize.getY());

    // Sensitive detector doesn't matter which logical volume is used
    G4VSensitiveDetector *sensDet=new CML2DummySD("Dummy ROG phantom");
    ROPhantomYDivisionLog->SetSensitiveDetector(sensDet);

    return ROPhyVol;
}
