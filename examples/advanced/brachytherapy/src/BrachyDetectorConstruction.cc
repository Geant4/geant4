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

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName):
  MaterialBox(0),DefaultMaterial(0),Iodio(0),Gold(0),titanium(0)
{
  NumVoxelX=300;
  NumVoxelZ=300;
  m_BoxDimX=30*cm ;
  m_BoxDimY=30*cm;
  m_BoxDimZ=30*cm;
  theUserLimitsForWater = NULL; 
  ComputeDimVoxel();
 
     

// set fUserLimit true to have the limit on step lenght
 //Geant mi da' l'energia depositata su uno step.
 fUseUserLimits = false;//non fa nulla se c'e false ,devo mettere true 
 theMaxStepInWater = .5*mm;
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
 /*

 A = 207.19*g/mole;
 d = 11.35*g/cm3;
 G4Material* matPb = new G4Material("Lead",Z=82.,A,d);
 
 */
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

 //gold(chimica degli elementi N.N Greenwood,A.Earnshaw)
 A=196.97*g/mole;
 d=19.32*g/cm3;
 G4Material* gold=new G4Material("gold",Z=79.,A,d);

 // Iridium (Medical Physics, Vol 25, No 10, Oct 1998)
 /*
 d = 22.42*g/cm3;
 A = 191.96260*g/mole ;
 G4Material* matIr192 = new G4Material("Iridium",Z=77.,A,d);
 */

 //IodiumCore(chimica degli elementi N.N Greenwood,A.Earnshaw)
 A=124.9*g/mole;
 d=4.862*g/cm3;
 G4Material* matI=new G4Material("Iodium",Z=53.,A,d);

 //ceramica(Medical Physics, May 2000)
  
  d=2.88*g/cm3;
 G4Material* ceramica=new G4Material("allumina",d,2);
 ceramica->AddElement(elAl,2);
 ceramica->AddElement(elO,3);

 // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
 d = 8.02*g/cm3 ;
 G4Material* matSteel = new G4Material("Stainless steel",d,5);
 matSteel->AddElement(elMn, 0.02);
 matSteel->AddElement(elSi, 0.01);
 matSteel->AddElement(elCr, 0.19);
 matSteel->AddElement(elNi, 0.10);
 matSteel->AddElement(elFe, 0.68);
 
 //G4cout << *(G4Material::GetMaterialTable()) << G4endl;



 MaterialBox=matH2O;
 DefaultMaterial=matAir;
 Iodio= matI;
 Gold=gold;
 titanium=Titanium;
}

G4VPhysicalVolume* BrachyDetectorConstruction::ConstructDetector()
{// Volumes
  ComputeDimVoxel();
 // EXPERIMENTAL HALL (our world volume)
 G4double ExpHall_x = 4.0*m;
 G4double ExpHall_y = 4.0*m;
 G4double ExpHall_z = 4.0*m;


  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  // Rotation Matrix
  G4RotationMatrix *rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateY(180.*deg);

  //World volume
 G4Box* ExpHall = new G4Box("ExpHall",ExpHall_x,ExpHall_y,ExpHall_z);
 G4LogicalVolume* ExpHallLog = new G4LogicalVolume(ExpHall,DefaultMaterial,"ExpHallLog",0,0,0);
 G4VPhysicalVolume* ExpHallPhys = new G4PVPlacement
			(0,G4ThreeVector(),"ExpHallPhys",ExpHallLog,0,false,0);

 // Water Box
 G4Box* WaterBox = new G4Box("WaterBox",m_BoxDimX/2,m_BoxDimY/2,m_BoxDimZ/2);
 WaterBoxLog = new G4LogicalVolume(WaterBox,MaterialBox,"WaterBoxLog",0,0,0);
 
 // create UserLimits
  if (theUserLimitsForWater != NULL) delete theUserLimitsForWater;
  theUserLimitsForWater = new G4UserLimits(//DBL_MAX,  //step max
					      //DBL_MAX,  // track max
					      theMaxStepInWater );
					     
					
 

 // attach UserLimits   
  if (fUseUserLimits) {
    WaterBoxLog->SetUserLimits(theUserLimitsForWater);
  }






 G4VPhysicalVolume* WaterBoxPhys = new G4PVPlacement
			(0,G4ThreeVector(),WaterBoxLog,"WaterBoxPhys",ExpHallLog,false,0);

 //Barbatrucco
 G4Tubs* DefaultTub=new G4Tubs("DefaultTub",
                               0.*mm,
                               0.40*mm,
                               1.84*mm,
                               0.*deg,
                               360.*deg);
 G4LogicalVolume* DefaultTub_log = new G4LogicalVolume(DefaultTub,DefaultMaterial,"DefaultTub_Log");
 G4VPhysicalVolume* DefaultTub_Phys = new G4PVPlacement
			(0,G4ThreeVector(0,0,0),DefaultTub_log,
			"DefaultTub_Log",WaterBoxLog,false,0); 
 //  Capsule main body
 G4double CapsuleR=0.35*mm;
 G4Tubs* Capsule = new G4Tubs("Capsule",
                              CapsuleR,//raggio interno
                              0.40*mm,//raggio esterno
                              1.84*mm,//meta' lunghezza 
                              0.*deg,//angolo di partenza
                              360.*deg);//angolo di rotazione
 
 G4LogicalVolume* CapsuleLog = new G4LogicalVolume(Capsule,titanium,"CapsuleLog");
 G4VPhysicalVolume* CapsulePhys = new G4PVPlacement
			(0,G4ThreeVector(0,0,0),CapsuleLog,
			"CapsulePhys",DefaultTub_log,false,0);
 
 // Capsule tip

 G4Sphere* CapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.40*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  G4LogicalVolume* CapsuleTipLog = new G4LogicalVolume(CapsuleTip,titanium,"CapsuleTipLog");
 G4VPhysicalVolume* CapsuleTipPhys1 = new G4PVPlacement
			(0,G4ThreeVector(0.,0.,1.84*mm),CapsuleTipLog,
			"CapsuleTipPhys",WaterBoxLog,false,0);

 G4VPhysicalVolume* CapsuleTipPhys2 = new G4PVPlacement
			(rotateMatrix,G4ThreeVector(0.,0.,-1.84*mm),CapsuleTipLog,"CapsuleTipPhys",WaterBoxLog,false,0);
 
// I-125 core
 // al centro del sistema di rif
  G4Tubs* IodiumCore = new G4Tubs("ICore",0.085*mm,0.30*mm,1.75*mm,0.*deg,360.*deg);
  
 G4LogicalVolume* ICoreLog = new G4LogicalVolume(IodiumCore,Iodio,"IridiumCoreLog");
 G4VPhysicalVolume* ICorePhys = new G4PVPlacement
			(0,G4ThreeVector(0.,0.,0.),ICoreLog,
			"IodiumCorePhys",DefaultTub_log,false,0);
 
 //gold marker al centro del sistema di rif
 

 G4Tubs* Marker=new G4Tubs("GoldenMarker",0.*mm,0.085*mm,1.75*mm,
                                                      0.*deg,360.*deg);
 G4LogicalVolume* MarkerLog = new G4LogicalVolume(Marker,Gold,"MarkerLog");
 G4VPhysicalVolume* MarkerPhys=new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),MarkerLog,"MarkerPhys",DefaultTub_log,false,0);
 
 // Sensitive Detector and ReadOut geometry definition
 //
 G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
 BrachyWaterBoxSD* pWaterBoxSD = new BrachyWaterBoxSD(m_SDName,NumVoxelX,NumVoxelZ);//this con GUi cambiare il WaterBox 
 if(pWaterBoxSD)
	{
	G4String ROGeometryName = "WaterBoxROGeometry";
	BrachyWaterBoxROGeometry* pWaterBoxROGeometry = new BrachyWaterBoxROGeometry(ROGeometryName,m_BoxDimX,m_BoxDimZ,NumVoxelX,NumVoxelZ);
	pWaterBoxROGeometry->BuildROGeometry();
	pWaterBoxSD->SetROgeometry(pWaterBoxROGeometry);
	pSDManager->AddNewDetector(pWaterBoxSD);
	
        WaterBoxLog->SetSensitiveDetector(pWaterBoxSD);
	//   CapsuleLog->SetSensitiveDetector(pWaterBoxSD);
        //CapsuleTipLog->SetSensitiveDetector(pWaterBoxSD);
	// DefaultLog->SetSensitiveDetector(pWaterBoxSD);
	 //    DefaultLog2->SetSensitiveDetector(pWaterBoxSD);
        //ICoreLog->SetSensitiveDetector(pWaterBoxSD);
        //MarkerLog->SetSensitiveDetector(pWaterBoxSD);
	}

   //                                        
  // Visualization attributes
 ExpHallLog->SetVisAttributes (G4VisAttributes::Invisible);
  
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(lblue);
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxVisAtt->SetForceWireframe(true);

  WaterBoxLog->SetVisAttributes(simpleBoxVisAtt);
  
   G4VisAttributes* simpleMarkerVisAtt= new G4VisAttributes(lblue);
 simpleMarkerVisAtt->SetVisibility(true);
  simpleMarkerVisAtt->SetForceSolid(true);
  MarkerLog->SetVisAttributes(simpleMarkerVisAtt);
  //  DefaultTub_log->SetVisAttributes(simpleMarkerVisAtt);
  
  
G4VisAttributes* simpleIodiumVisAtt= new G4VisAttributes(magenta);
    simpleIodiumVisAtt->SetVisibility(true);
  simpleIodiumVisAtt->SetForceWireframe(true);
  ICoreLog->SetVisAttributes(simpleIodiumVisAtt);
   
  
  
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
        <<"material of the box : "
        <<MaterialBox->GetName()
        <<G4endl
	<<"the source is at the centre of the detector"
         <<G4endl
	 

<<"-------------------------------------------------------------------------"
	 << G4endl;
}

//da sistemare
void BrachyDetectorConstruction::SetMaterial(G4String materialChoice)
{
 G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {MaterialBox = pttoMaterial;
      WaterBoxLog->SetMaterial(pttoMaterial); 
      PrintDetectorParameters();

     }
}

void BrachyDetectorConstruction::SetDimension(G4double val)

{ 
 m_BoxDimX=val ;
 m_BoxDimY=val ;
 m_BoxDimZ=val ;

}


void BrachyDetectorConstruction::SetNumVoxel(G4int val)
{ 
  NumVoxelX=val;
  NumVoxelZ=val;
  
       
}
void  BrachyDetectorConstruction::SetMaxStepInWater(G4double value)
{ 
  theMaxStepInWater = value; 
  if (theUserLimitsForWater != NULL) 
  {
    theUserLimitsForWater->SetMaxAllowedStep(value);
  }
}

void  BrachyDetectorConstruction::UseUserLimits(G4bool isUse)
{
  fUseUserLimits = isUse;
  if ( fUseUserLimits && (theUserLimitsForWater!= NULL)) 
  {WaterBoxLog->SetUserLimits(theUserLimitsForWater);
  }    
} 

void BrachyDetectorConstruction::UpdateGeometry()
{ 
   PrintDetectorParameters();
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructDetector());
}













