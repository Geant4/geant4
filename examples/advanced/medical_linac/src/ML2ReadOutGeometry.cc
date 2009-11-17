#include "ML2ReadOutGeometry.h"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "ML2DummySD.h"
#include "G4PVReplica.hh"

CML2ReadOutGeometry::CML2ReadOutGeometry(const G4RotationMatrix *m, G4ThreeVector *v):ROPhyVol(0)
{
	// Build the world volume 
	G4RotationMatrix *ml=new G4RotationMatrix;
	*ml=*m;
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4ThreeVector halfSizeWorld, centre;
	centre.set(0.*mm, 0.*mm, 0.*mm); centre+=*v;
	halfSizeWorld.set(3000.*mm, 3000*mm, 3000*mm);
	G4Box *ROphmWorldB = new G4Box("ROphmWorldG", halfSizeWorld.getX(), halfSizeWorld.getY(), halfSizeWorld.getZ());
	G4LogicalVolume *ROphmWorldLV = new G4LogicalVolume(ROphmWorldB, Vacuum, "ROphmWorldL", 0, 0, 0);
	this->ROPhyVol= new G4PVPlacement(ml, centre, "ROphmWorldPV", ROphmWorldLV, 0, false, 0);
}

CML2ReadOutGeometry::~CML2ReadOutGeometry(void)
{
	delete this->ROPhyVol;
}
void CML2ReadOutGeometry::setBuildData(G4ThreeVector centre, G4ThreeVector halfSize, G4int NumberOfVoxelsAlongX, G4int NumberOfVoxelsAlongY, G4int NumberOfVoxelsAlongZ)
{
	this->centre=centre;
	this->halfSize=halfSize;
	this->NumberOfVoxelsAlongX=NumberOfVoxelsAlongX;
	this->NumberOfVoxelsAlongY=NumberOfVoxelsAlongY;
	this->NumberOfVoxelsAlongZ=NumberOfVoxelsAlongZ;
}
G4VPhysicalVolume* CML2ReadOutGeometry::Build()
{
	// Build RO Zone
	G4Material *Vacuum=G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Box *ROBox = new G4Box("ROBox", this->halfSize.getX(), this->halfSize.getY(), this->halfSize.getZ());
	G4LogicalVolume *ROLV = new G4LogicalVolume(ROBox, Vacuum, "ROLV", 0, 0, 0);
	G4VPhysicalVolume *ROPV = new G4PVPlacement(0, this->centre, "ROPV", ROLV, this->ROPhyVol, false, 0);

 // ROGeomtry: Voxel division
 
	G4double halfXVoxelDimensionX, halfXVoxelDimensionY, halfXVoxelDimensionZ;

	halfXVoxelDimensionX=this->halfSize.getX()/this->NumberOfVoxelsAlongX;
	halfXVoxelDimensionY=this->halfSize.getY()/this->NumberOfVoxelsAlongY;
	halfXVoxelDimensionZ=this->halfSize.getZ()/this->NumberOfVoxelsAlongZ;

	G4double voxelXThicknessX = 2*halfXVoxelDimensionX;
	G4double voxelXThicknessY = 2*halfXVoxelDimensionY;
	G4double voxelXThicknessZ = 2*halfXVoxelDimensionZ;


  // X division first... slice along X axis
  G4Box *ROPhantomXDivision = new G4Box("ROPhantomXDivision",
					halfXVoxelDimensionX,
					this->halfSize.getY(),
					this->halfSize.getZ());

  G4LogicalVolume *ROPhantomXDivisionLog = new G4LogicalVolume(ROPhantomXDivision,
							       Vacuum,
							       "ROPhantomXDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomXDivisionPhys = new G4PVReplica("ROPhantomXDivisionPhys",
                                                              ROPhantomXDivisionLog,
								ROPV,
                                                              kXAxis,
                                                              NumberOfVoxelsAlongX,
                                                              voxelXThicknessX,
								-this->halfSize.getX());
  // ...then Z division
  
  G4Box *ROPhantomZDivision = new G4Box("ROPhantomZDivision",
					halfXVoxelDimensionX,
					this->halfSize.getY(), 
					halfXVoxelDimensionZ);

  G4LogicalVolume *ROPhantomZDivisionLog = new G4LogicalVolume(ROPhantomZDivision,
							       Vacuum,
							       "ROPhantomZDivisionLog",
							       0,0,0);

  G4VPhysicalVolume *ROPhantomZDivisionPhys = new G4PVReplica("ROPhantomZDivisionPhys",
							      ROPhantomZDivisionLog,
							      ROPhantomXDivisionPhys,
							      kZAxis,
							      this->NumberOfVoxelsAlongZ,
							      voxelXThicknessZ,
								  -this->halfSize.getZ());
  // ...then Y  division

  G4Box *ROPhantomYDivision = new G4Box("ROPhantomYDivision",
					halfXVoxelDimensionX, 
					halfXVoxelDimensionY,
					halfXVoxelDimensionZ);

  G4LogicalVolume *ROPhantomYDivisionLog = new G4LogicalVolume(ROPhantomYDivision,
							       Vacuum,
							       "ROPhantomYDivisionLog",
							       0,0,0);
 G4VPhysicalVolume *ROPhantomYDivisionPhys = 0;
 ROPhantomYDivisionPhys = new G4PVReplica("ROPhantomYDivisionPhys",
							      ROPhantomYDivisionLog,
							      ROPhantomZDivisionPhys,
							      kYAxis,
							      this->NumberOfVoxelsAlongY,
							      voxelXThicknessY,
								  -this->halfSize.getY());

	// Sensitive detector doesn't matter which logical volume is used 
	G4VSensitiveDetector *sensDet=new CML2DummySD("Dummy ROG phantom");
	ROPhantomYDivisionLog->SetSensitiveDetector(sensDet);
  return ROPhyVol;
}
