//    ************************************
//    *                                  *
//    *   BrachyWaterBoxROGeometry.cc    *
//    *                                  *
//    ************************************

#include "BrachyWaterBoxROGeometry.hh"
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

//....

BrachyWaterBoxROGeometry::BrachyWaterBoxROGeometry(G4String aString,G4double DetDimX,G4double DetDimZ,G4int NumVoxelX,G4int NumVoxelZ)
: G4VReadOutGeometry(aString),m_DetDimX(DetDimX),m_DetDimZ(DetDimZ),m_NumVoxelX(NumVoxelX),m_NumVoxelZ(NumVoxelZ)
{
}

//....

BrachyWaterBoxROGeometry::~BrachyWaterBoxROGeometry()
{
}

//....

G4VPhysicalVolume* BrachyWaterBoxROGeometry::Build()
{
 // A dummy material is used to fill the volumes of the readout geometry.
 // (It will be allowed to set a NULL pointer in volumes of such virtual
 // division in future, since this material is irrelevant for tracking.)

 G4Material* dummyMat  = new G4Material(name="dummyMat", 1., 1.*g/mole, 1.*g/cm3);

 // Slice thickness is the average of Voxel X and Z sizes
 G4double DetVoxel_y = (m_DetDimX/m_NumVoxelX+m_DetDimZ/m_NumVoxelZ)/2.0;

 G4double ExpHall_x = 4.0*m;
 G4double ExpHall_y = 4.0*m;
 G4double ExpHall_z = 4.0*m;

 G4double Det_x = m_DetDimX/2;
 G4double Det_y = DetVoxel_y;
 G4double Det_z = m_DetDimZ/2;

 G4double DetVoxelX_x = Det_x/m_NumVoxelX;
 G4double DetVoxelX_y = DetVoxel_y; 
 G4double DetVoxelX_z = Det_z;
 G4double DetVoxelX_dx = 2*DetVoxelX_x;

 G4double DetVoxelZ_x = Det_x;
 G4double DetVoxelZ_y = DetVoxel_y; 
 G4double DetVoxelZ_z = Det_z/m_NumVoxelZ;
 G4double DetVoxelZ_dz = 2*DetVoxelZ_z;

 G4Box *ROExpHall = new G4Box("ROExpHall",ExpHall_x,ExpHall_y,ExpHall_z);
 G4LogicalVolume *ROExpHallLog = new G4LogicalVolume(ROExpHall,dummyMat,"ROExpHallLog",0,0,0);
 G4VPhysicalVolume *ROExpHallPhys = new G4PVPlacement(0,G4ThreeVector(),"ROExpHallPhys",ROExpHallLog,0,false,0);
  
 G4Box *RODetector = new G4Box("RODetector", Det_x, Det_y, Det_z);
 G4LogicalVolume *RODetectorLog = new G4LogicalVolume(RODetector,dummyMat,"RODetectorLog",0,0,0);
 G4VPhysicalVolume *RODetectorPhys = new G4PVPlacement(0,G4ThreeVector(),"DetectorPhys",RODetectorLog,ROExpHallPhys,false,0);

 // ReadOut Voxel division
 
 // X division first...
  
 G4Box *RODetectorXDivision = new G4Box("RODetectorXDivision",DetVoxelX_x,DetVoxelX_y,DetVoxelX_z);
 G4LogicalVolume *RODetectorXDivisionLog = new G4LogicalVolume(RODetectorXDivision,dummyMat,"RODetectorXDivisionLog",0,0,0);
 G4VPhysicalVolume *RODetectorXDivisionPhys = new G4PVReplica("RODetectorXDivisionPhys",RODetectorXDivisionLog,RODetectorPhys,kXAxis,m_NumVoxelX,DetVoxelX_dx);

// ...then Z division
  
 G4Box *RODetectorZDivision = new G4Box("RODetectorZDivision",DetVoxelZ_x,DetVoxelZ_y,DetVoxelZ_z);
 G4LogicalVolume *RODetectorZDivisionLog = new G4LogicalVolume(RODetectorZDivision,dummyMat,"RODetectorZDivisionLog",0,0,0);
 G4VPhysicalVolume *RODetectorZDivisionPhys = new G4PVReplica("RODetectorZDivisionPhys",RODetectorZDivisionLog,RODetectorXDivisionPhys,kZAxis,m_NumVoxelZ,DetVoxelZ_dz);
  
 BrachyDummySD *dummySD = new BrachyDummySD;
 RODetectorZDivisionLog->SetSensitiveDetector(dummySD);

 return ROExpHallPhys;
}



