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
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName):
  NumVoxelX(0),NumVoxelZ(0),m_BoxDimX(0), m_BoxDimY(0), m_BoxDimZ(0),
  ExpHall(0),ExpHallLog(0),ExpHallPhys(0),
  WaterBox(0), WaterBoxLog(0), WaterBoxPhys(0),
   Capsule(0),CapsuleLog(0),CapsulePhys(0),
  CapsuleTip(0),CapsuleTipLog(0),CapsuleTipPhys(0),
  IridiumCore(0),IridiumCoreLog(0),IridiumCorePhys(0)
 
{
  NumVoxelX=300;
  NumVoxelZ=300;

  m_BoxDimX=30*cm ;
  m_BoxDimY=30*cm;
  m_BoxDimZ=30*cm;

  ComputeDimVoxel();
 
     
  m_SDName = SDName;//pointer to sensitive detector
 
}


//....

BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 
 
}

//....

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{ 
  DefineMaterials();
  return ConstructDetector();
}

void BrachyDetectorConstruction::DefineMaterials()
{
 // Define required materials

 G4double A;  // atomic mass
 G4double Z;  // atomic number
 G4double d;  // density
 
 // General elements
 
 A = 1.01*g/mole;
 G4Element* elH = new G4Element ("Hydrogen","H",Z=1.,A);
  
 A = 14.01*g/mole;
 G4Element* elN = new G4Element("Nitrogen","N",Z=7.,A);

 A = 16.00*g/mole;
  G4Element* elO = new G4Element("Oxygen","O",Z=8.,A);
   
 // Elements for source capsule and cable
 
  A=26.98*g/mole;
  G4Element* elAl=new G4Element("Aluminum","Al", Z=13.,A);


 A = 54.94*g/mole;
 G4Element* elMn  = new G4Element("Manganese","Mn",Z=25.,A);
 
 A = 28.09*g/mole;
 G4Element* elSi  = new G4Element("Silicon","Si",Z=14.,A);

 A = 52.00*g/mole;
 G4Element* elCr  = new G4Element("Chromium","Cr",Z=24.,A);

 A = 58.70*g/mole;
 G4Element* elNi  = new G4Element("Nickel","Ni",Z=28.,A);

 A = 55.85*g/mole;
 G4Element* elFe  = new G4Element("Iron","Fe",Z=26.,A);

// Lead material
 A = 207.19*g/mole;
 Z = 82;
 d = 11.35*g/cm3;
 G4Material* matPb = new G4Material("Lead",Z,A,d);

 // Iridium (Medical Physics, Vol 25, No 10, Oct 1998)
 d = 22.42*g/cm3;
 A = 191.96260*g/mole ;
 Z = 77;
 G4Material* matIr192 = new G4Material("Iridium",Z,A,d);

 //titanium
 A=47.88*g/mole;
 d=4.50*g/cm3;
  G4Material* Titanium=new G4Material("titanium" ,Z=22.,A,d);
 
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

 // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
 d = 8.02*g/cm3 ;
 G4Material* matSteel = new G4Material("Stainless steel",d,5);
 matSteel->AddElement(elMn, 0.02);
 matSteel->AddElement(elSi, 0.01);
 matSteel->AddElement(elCr, 0.19);
 matSteel->AddElement(elNi, 0.10);
 matSteel->AddElement(elFe, 0.68);
 
//--elements for Iodium source which will be introduced in the next release

 //gold(chimica degli elementi N.N Greenwood,A.Earnshaw)
 A=196.97*g/mole;
 d=19.32*g/cm3;
 G4Material* gold=new G4Material("gold",Z=79.,A,d);


 //IodiumCore(chimica degli elementi N.N Greenwood,A.Earnshaw)
 A=124.9*g/mole;
 d=4.862*g/cm3;
 G4Material* matI=new G4Material("Iodium",Z=53.,A,d);

 //ceramica(Medical Physics, May 2000)
  
  d=2.88*g/cm3;
 G4Material* ceramica=new G4Material("allumina",d,2);
 ceramica->AddElement(elAl,2);
 ceramica->AddElement(elO,3);


 water=matH2O;
 air= matAir;
 CapsuleMat=matSteel;
 IridiumMat=matIr192;

}

G4VPhysicalVolume* BrachyDetectorConstruction::ConstructDetector()
{// Volumes
  ComputeDimVoxel();
 // EXPERIMENTAL HALL (our world volume)
 G4double ExpHall_x = 4.0*m;
 G4double ExpHall_y = 4.0*m;
 G4double ExpHall_z = 4.0*m;

  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.75, .75, .75) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75);

 
 //World volume
  ExpHall = new G4Box("ExpHall",ExpHall_x,ExpHall_y,ExpHall_z);
  ExpHallLog = new G4LogicalVolume(ExpHall,air,"ExpHallLog",0,0,0);
 ExpHallPhys = new G4PVPlacement(0,G4ThreeVector(),"ExpHallPhys",ExpHallLog,NULL,false,0);

 // Water Box
 WaterBox = new G4Box("WaterBox",m_BoxDimX/2,m_BoxDimY/2,m_BoxDimZ/2);
 WaterBoxLog = new G4LogicalVolume(WaterBox,water,"WaterBoxLog",0,0,0);
 


WaterBoxPhys = new G4PVPlacement(0,G4ThreeVector(),"WaterBoxPhys",WaterBoxLog,ExpHallPhys,false,0);



// Capsule main body

 G4Tubs* Capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
 G4LogicalVolume* CapsuleLog = new G4LogicalVolume(Capsule,CapsuleMat,"CapsuleLog");
 G4VPhysicalVolume* CapsulePhys = new G4PVPlacement(0,G4ThreeVector(0,0,-1.975),"CapsulePhys",CapsuleLog,WaterBoxPhys,false,0);

 // Capsule tip

 G4Sphere* CapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
 G4LogicalVolume* CapsuleTipLog = new G4LogicalVolume(CapsuleTip,CapsuleMat,"CapsuleTipLog");
 G4VPhysicalVolume* CapsuleTipPhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,1.75*mm),"CapsuleTipPhys",CapsuleTipLog,WaterBoxPhys,false,0);

 // Iridium core

 G4Tubs* IridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
 G4LogicalVolume* IridiumCoreLog = new G4LogicalVolume(IridiumCore,IridiumMat,"IridiumCoreLog");
 G4VPhysicalVolume* IridiumCorePhys = new G4PVPlacement(0,G4ThreeVector(),"IridiumCorePhys",IridiumCoreLog,CapsulePhys,false,0);
	


/*

  DefaultTub=new G4Tubs("DefaultTub",
                               0.*mm,
                               0.40*mm,
                               1.84*mm,
                               0.*deg,
                               360.*deg);
 DefaultTub_log = new G4LogicalVolume(DefaultTub,DefaultMaterial,"DefaultTub_Log");
 DefaultTub_Phys = new G4PVPlacement(0,G4ThreeVector(),"DefaultTub_Phys",DefaultTub_log,
			 WaterBoxPhys,false,0); 
 //  Capsule main body 

// Rotation Matrix
  G4RotationMatrix*  rotateMatrix=new G4RotationMatrix(); //=new G4RotationMatrix() ;
  rotateMatrix->rotateX(180.0*deg);
 G4double CapsuleR=0.35*mm;
  Capsule = new G4Tubs("Capsule",
                              CapsuleR,//raggio interno
                              0.40*mm,//raggio esterno
                              1.84*mm,//meta' lunghezza 
                              0.*deg,//angolo di partenza
                              360.*deg);//angolo di rotazione
 
 CapsuleLog = new G4LogicalVolume(Capsule,titanium,"CapsuleLog");
 CapsulePhys = new G4PVPlacement
			(0,G4ThreeVector(),"CapsulePhys",CapsuleLog,
			DefaultTub_Phys,false,0);
 
 // Capsule tip

  CapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.40*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  CapsuleTipLog = new G4LogicalVolume(CapsuleTip,titanium,"CapsuleTipLog");
 CapsuleTipPhys1 = new G4PVPlacement
	      (0,G4ThreeVector(0.,0.,1.84*mm), "CapsuleTipPhys1",CapsuleTipLog,WaterBoxPhys,false,0);


 CapsuleTipPhys2 = new G4PVPlacement
			( rotateMatrix, G4ThreeVector(0,0,-1.84*mm),"CapsuleTipPhys2",CapsuleTipLog,WaterBoxPhys,false,0);
 
  CapsuleTipPhys2 = new G4PVPlacement
			( G4Transform3D(rotateMatrix, G4ThreeVector(0,0,-1.84*mm)),"CapsuleTipPhys2",CapsuleTipLog,WaterBoxPhys,false,0);
 
// I-125 core

   IodiumCore = new G4Tubs("ICore",0.085*mm,0.35*mm,1.75*mm,0.*deg,360.*deg);
  
    ICoreLog = new G4LogicalVolume(IodiumCore,Iodio,"IodiumCoreLog");
    ICorePhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
			"IodiumCorePhys",ICoreLog,DefaultTub_Phys,false,0);
 
 //gold marker al centro del sistema di rif
 

  Marker=new G4Tubs("GoldenMarker",0.*mm,0.085*mm,1.75*mm,
                                                      0.*deg,360.*deg);
  MarkerLog = new G4LogicalVolume(Marker,Gold,"MarkerLog");
  MarkerPhys=new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),"MarkerPhys",MarkerLog,DefaultTub_Phys,false,0);
 */
 // Sensitive Detector and ReadOut geometry definition
 //
 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 BrachyWaterBoxSD* pWaterBoxSD = new BrachyWaterBoxSD(m_SDName,NumVoxelX,NumVoxelZ);
 if(pWaterBoxSD)
	{
	G4String ROGeometryName = "WaterBoxROGeometry";
	BrachyWaterBoxROGeometry* pWaterBoxROGeometry = new BrachyWaterBoxROGeometry(ROGeometryName,m_BoxDimX,m_BoxDimZ,NumVoxelX,NumVoxelZ);
	pWaterBoxROGeometry->BuildROGeometry();
	pWaterBoxSD->SetROgeometry(pWaterBoxROGeometry);
	pSDManager->AddNewDetector(pWaterBoxSD);
	
         WaterBoxLog->SetSensitiveDetector(pWaterBoxSD);
	  CapsuleLog->SetSensitiveDetector(pWaterBoxSD);
          CapsuleTipLog->SetSensitiveDetector(pWaterBoxSD);
	  IridiumCoreLog->SetSensitiveDetector(pWaterBoxSD);
        
	}

   //                                        
  // Visualization attributes
 //
 ExpHallLog->SetVisAttributes (G4VisAttributes::Invisible);
  
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(lblue);
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxVisAtt->SetForceWireframe(true);

  WaterBoxLog->SetVisAttributes(simpleBoxVisAtt);
  
   
  
  G4VisAttributes* simpleIridiumVisAtt= new G4VisAttributes(magenta);
    simpleIridiumVisAtt->SetVisibility(true);
    simpleIridiumVisAtt->SetForceWireframe(true);
  IridiumCoreLog->SetVisAttributes(simpleIridiumVisAtt);
   
  
  
 G4VisAttributes* simpleCapsuleVisAtt= new G4VisAttributes(red);
 simpleCapsuleVisAtt->SetVisibility(true);  
 simpleCapsuleVisAtt->SetForceWireframe(true);
 CapsuleLog->SetVisAttributes( simpleCapsuleVisAtt);
  
   G4VisAttributes* simpleCapsuleTipVisAtt= new G4VisAttributes(red);
 simpleCapsuleTipVisAtt->SetVisibility(true);  
 simpleCapsuleTipVisAtt->SetForceSolid(true);
 CapsuleTipLog->SetVisAttributes( simpleCapsuleTipVisAtt);
  
 
 //always return the physical World
  
 PrintDetectorParameters();
 return ExpHallPhys;
}


void BrachyDetectorConstruction::PrintDetectorParameters()
{
 G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	<<"the detector is a  box whose size is: " 
        <<G4endl
        <<m_BoxDimX/cm
        << " cm * "
        <<m_BoxDimY/cm
        << " cm * "
        <<m_BoxDimZ/cm
        << " cm"
        << G4endl
        <<"numVoxel: "
        <<NumVoxelX
	<<G4endl
	<<"dim voxel: "
        <<dimVoxel/mm
        <<"mm"
	<<G4endl 
        <<"material of the box :water "
        <<G4endl
	<<"the source is at the center  of the detector"
         <<G4endl
	 

<<"-------------------------------------------------------------------------"
	 << G4endl;
}
















