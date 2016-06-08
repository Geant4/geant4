//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.cc     *
//    *                                      *
//    ****************************************

#include "BrachyWaterBoxROGeometry.hh"
#include "BrachyWaterBoxSD.hh"
#include "BrachyDetectorConstruction.hh"

#include "G4CSGSolid.hh"
#include "G4Sphere.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName,G4int NumVoxelX,G4int NumVoxelZ) :
	m_NumVoxelX(NumVoxelX),m_NumVoxelZ(NumVoxelZ),m_BoxDimX(30*cm),m_BoxDimY(30*cm),m_BoxDimZ(30*cm)
{
 m_SDName = SDName;
}

//....

BrachyDetectorConstruction::~BrachyDetectorConstruction()
{
}

//....

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{
 // Define required materials

 G4double A;  // atomic mass
 G4double Z;  // atomic number
 G4double d;  // density
 
 // General elements
 
 A = 1.01*g/mole;
 Z = 1;
 G4Element* elH = new G4Element ("Hydrogen","H",Z,A);
  
 A = 14.01*g/mole;
 Z = 7;
 G4Element* elN = new G4Element("Nitrogen","N",Z,A);

 A = 16.00*g/mole;
 Z = 8;
 G4Element* elO = new G4Element("Oxygen","O",Z,A);
   
 // Elements for source capsule and cable

 A = 54.94*g/mole;
 Z = 25;
 G4Element* elMn  = new G4Element("Manganese","Mn",Z,A);
 
 A = 28.09*g/mole;
 Z = 14;
 G4Element* elSi  = new G4Element("Silicon","Si",Z,A);

 A = 52.00*g/mole;
 Z = 24;
 G4Element* elCr  = new G4Element("Chromium","Cr",Z,A);

 A = 58.70*g/mole;
 Z = 28;
 G4Element* elNi  = new G4Element("Nickel","Ni",Z,A);

 A = 55.85*g/mole;
 Z = 26;
 G4Element* elFe  = new G4Element("Iron","Fe",Z,A);

 // Lead material
 A = 207.19*g/mole;
 Z = 82;
 d = 11.35*g/cm3;
 G4Material* matPb = new G4Material("Lead",Z,A,d);

 // Air material
 d = 1.290*mg/cm3;
 G4Material* matAir = new G4Material("Air",d,2);
 matAir->AddElement(elN,0.7);
 matAir->AddElement(elO,0.3);

 // Water
 d = 1.000*g/cm3;
 G4Material* matH2O = new G4Material("Water",d,2);
 matH2O->AddElement(elH,2);
 matH2O->AddElement(elO,1);

 // Iridium (Medical Physics, Vol 25, No 10, Oct 1998)
 d = 22.42*g/cm3;
 A = 191.96260*g/mole ;
 Z = 77;
 G4Material* matIr192 = new G4Material("Iridium",Z,A,d);

 // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
 d = 8.02*g/cm3 ;
 G4Material* matSteel = new G4Material("Stainless steel",d,5);
 matSteel->AddElement(elMn, 0.02);
 matSteel->AddElement(elSi, 0.01);
 matSteel->AddElement(elCr, 0.19);
 matSteel->AddElement(elNi, 0.10);
 matSteel->AddElement(elFe, 0.68);
 
 // Volumes

 // EXPERIMENTAL HALL (our world volume)
 G4double ExpHall_x = 4.0*m;
 G4double ExpHall_y = 4.0*m;
 G4double ExpHall_z = 4.0*m;
  
 G4Box* ExpHall = new G4Box("ExpHall",ExpHall_x,ExpHall_y,ExpHall_z);
 G4LogicalVolume* ExpHallLog = new G4LogicalVolume(ExpHall,matAir,"ExpHallLog",0,0,0);
 G4VPhysicalVolume* ExpHallPhys = new G4PVPlacement(0,G4ThreeVector(),"ExpHallPhys",ExpHallLog,0,false,0);

 // Water Box

 G4Box* WaterBox = new G4Box("WaterBox",m_BoxDimX/2,m_BoxDimY/2,m_BoxDimZ/2);
 G4LogicalVolume* WaterBoxLog = new G4LogicalVolume(WaterBox,matH2O,"WaterBoxLog",0,0,0);
 G4VPhysicalVolume* WaterBoxPhys = new G4PVPlacement(0,G4ThreeVector(),WaterBoxLog,"WaterBoxPhys",ExpHallLog,false,0);
 
 // Capsule main body

 G4Tubs* Capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
 G4LogicalVolume* CapsuleLog = new G4LogicalVolume(Capsule,matSteel,"CapsuleLog");
 G4VPhysicalVolume* CapsulePhys = new G4PVPlacement(0,G4ThreeVector(0,0,-1.975),CapsuleLog,"CapsulePhys",WaterBoxLog,false,0);

 // Capsule tip

 G4Sphere* CapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
 G4LogicalVolume* CapsuleTipLog = new G4LogicalVolume(CapsuleTip,matSteel,"CapsuleTipLog");
 G4VPhysicalVolume* CapsuleTipPhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,1.75*mm),CapsuleTipLog,"CapsuleTipPhys",WaterBoxLog,false,0);

 // Iridium core

 G4Tubs* IridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
 G4LogicalVolume* IridiumCoreLog = new G4LogicalVolume(IridiumCore,matIr192,"IridiumCoreLog");
 G4VPhysicalVolume* IridiumCorePhys = new G4PVPlacement(0,G4ThreeVector(),IridiumCoreLog,"IridiumCorePhys",CapsuleLog,false,0);
	
 // Sensitive Detector and ReadOut geometry definition

 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 BrachyWaterBoxSD* pWaterBoxSD = new BrachyWaterBoxSD(m_SDName,m_NumVoxelX,m_NumVoxelZ);
 if(pWaterBoxSD)
	{
	G4String ROGeometryName = "WaterBoxROGeometry";
	BrachyWaterBoxROGeometry* pWaterBoxROGeometry = new BrachyWaterBoxROGeometry(ROGeometryName,m_BoxDimX,m_BoxDimZ,m_NumVoxelX,m_NumVoxelZ);
	pWaterBoxROGeometry->BuildROGeometry();
	pWaterBoxSD->SetROgeometry(pWaterBoxROGeometry);
	pSDManager->AddNewDetector(pWaterBoxSD);
	WaterBoxLog->SetSensitiveDetector(pWaterBoxSD);
	CapsuleLog->SetSensitiveDetector(pWaterBoxSD);
	CapsuleTipLog->SetSensitiveDetector(pWaterBoxSD);
	IridiumCoreLog->SetSensitiveDetector(pWaterBoxSD);
	}

 return ExpHallPhys;
}
