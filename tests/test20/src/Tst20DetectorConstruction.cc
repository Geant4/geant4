// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20DetectorConstruction.cc,v 1.1 2001-05-24 19:49:30 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ Tst20DetectorConstruction  ------
//
// ************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst20DetectorConstruction.hh"
#include "Tst20DetectorMessenger.hh"

#include "Tst20AnticoincidenceSD.hh"
#include "Tst20TrackerSD.hh"
#include "Tst20TrackerROGeometry.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20DetectorConstruction::Tst20DetectorConstruction()
  :solidWorld(0),logicWorld(0),physiWorld(0),
   solidPayload(0),logicPayload(0),physiPayload(0),
   solidTKR(0),logicTKR(0),physiTKR(0),
   solidACT(0),logicACT(0),physiACT(0),
   solidACL1(0),logicACL1(0),physiACL1(0),
   solidACL2(0),logicACL2(0),physiACL2(0),
   solidConverter(0),logicConverter(0),physiConverter(0),
   solidTKRDetector(0),logicTKRDetector(0),physiTKRDetector(0),
   solidPlane(0),logicPlane(0),physiPlane(0)
  
{
  // default parameter values of the payload
  
  ConverterThickness = 1.*micrometer;
  TKRSiliconThickness = 400.*micrometer;
  TKRSiliconTileXY = 9.*cm;
  TKRSiliconPitch = 1.*cm; 
  TKRLayerDistance = 3.*cm;
  TKRTileDistance = 1.*cm;
  SiliconGuardRing = 1.5*mm;
  TKRViewsDistance = 1.*mm;
  NbOfTKRLayers = 10;
  NbOfTKRTiles = 2;
  NbOfACDTopTiles = 2;
  NbOfACDLateralTiles = 4;
  ACDThickness = 1.*cm;
  ACDTKRDistance = 5.*cm;
  ComputePayloadParameters();

  // create commands for interactive definition of the payload
  detectorMessenger = new Tst20DetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20DetectorConstruction::~Tst20DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Tst20DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructPayload();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::DefineMaterials()
{ 
  
  G4String name, symbol;    
  G4double a, z, density;            
  
  G4int ncomponents, natoms;
  G4double abundance, fractionmass;
  G4double temperature, pressure;

  //
  // define Elements
  //
  
  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  
  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.006*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen"  ,symbol="N" , z= 7., a);
  
  a = 15.99*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 126.904*g/mole;
  G4Element* I  = new G4Element(name="Iodine"  ,symbol="I" , z= 53., a);

  a = 132.905*g/mole;
  G4Element* Cs  = new G4Element(name="Cesium"  ,symbol="Cs" , z= 55., a);

//
// define simple materials
//

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  density = 2.333*g/cm3;
  a = 28.09*g/mole;
  G4Material* Si = new G4Material(name="Silicon",z=14., a,density);
  
  density = 19.3*g/cm3;
  a = 183.84*g/mole;
  G4Material* W = new G4Material(name="Tungsten", z=74., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead", z=82., a, density);
  
  density = 7.87*g/cm3;
  a= 55.845*g/mole;
  G4Material* Fe = new G4Material(name="Iron", z=26.,a,density);
  
  //
  // define a material from elements.   case 1: chemical molecule
  //
  
  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  
  density = 4.53*g/cm3;
  G4Material* CsI = new G4Material(name="CesiumIodide", density, ncomponents=2);
  CsI->AddElement(C, natoms=5);
  CsI->AddElement(H, natoms=5);
  
  //
  // define a material from elements.   case 2: mixture by fractional mass
  //
  
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
//
// examples of vacuum
//

  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,kStateGas,temperature,pressure);
  
  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;         //from PhysicalConstants.h
  G4Material* beam = new G4Material(name="Beam", density, ncomponents=1,
				    kStateGas,temperature,pressure);
  beam->AddMaterial(Air, fractionmass=1.);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //default materials of the payload

  ConverterMaterial = Pb;
  defaultMaterial  = vacuum;
  ACDTMaterial = Sci;
  ACDMaterial = CsI;
  TKRMaterial = Si;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Tst20DetectorConstruction::ConstructPayload()
{
  // complete the Payload parameters definition 
  ComputePayloadParameters();
  
  //     
  // World
  //

  solidWorld = new G4Box("World",				     
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  
  
  logicWorld = new G4LogicalVolume(solidWorld,		
                                   defaultMaterial,	
                                   "World");		
  
  physiWorld = new G4PVPlacement(0,G4ThreeVector(),"World",logicWorld,
                                 0,false,0);			
 

  //                               
  // Payload
  //  
  
  solidPayload=0; logicPayload=0; physiPayload=0;
  solidTKR=0;logicTKR=0;physiTKR=0;
  solidACT=0;logicACT=0;physiACT=0;
  solidACL1=0;logicACL1=0;physiACL1=0; 
  solidACL2=0;logicACL2=0;physiACL2=0; 
  solidConverter=0;logicConverter=0;physiConverter=0; 
  solidTKRDetector=0;logicTKRDetector=0;physiTKRDetector=0;
  solidPlane=0;logicPlane=0;physiPlane=0;
  


  // Payload
  
  solidPayload = new G4Box("Payload",		
			   PayloadSizeXY/2,
			   PayloadSizeXY/2,
			   PayloadSizeZ/2);
  
  logicPayload = new G4LogicalVolume(solidPayload,	
				     defaultMaterial,	
				     "Payload");  
  
  physiPayload = new G4PVPlacement(0,		
				   G4ThreeVector(),	
				   "Payload",	
				   logicPayload,	
				   physiWorld,	
				   false,		
				   0);	       

  //                               
  // Anticoincidence Top (ACT)
  //
  
  solidACT = new G4Box("ACT",			
		       ACTSizeXY/2,ACTSizeXY/2,ACTSizeZ/2); 
  
  logicACT = new G4LogicalVolume(solidACT,ACDTMaterial,"ACT");
      
  physiACT = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+TKRSizeZ
                                             +2*ACDTKRDistance+ACTSizeZ+
                                             ACTSizeZ/2),
			       "ACT",		
			       logicACT,
			       physiPayload,
			       false,	
			       0);	
    
  //                                 
  // Anticoincidence Bottom (ACB)
  //
  
  solidACB = new G4Box("ACB",			
		       ACTSizeXY/2,ACTSizeXY/2,ACTSizeZ/2); 
  
  logicACB = new G4LogicalVolume(solidACB,ACDMaterial,"ACB");
      
  physiACB = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+ACTSizeZ/2),
			       "ACB",		
			       logicACB,
			       physiPayload,
			       false,	
			       0);	
  
  //                                 
  // Payload (TKR)
  //
  
  solidTKR = new G4Box("TKR",			
		       TKRSizeXY/2,TKRSizeXY/2,TKRSizeZ/2); 
  
  logicTKR = new G4LogicalVolume(solidTKR,
				 defaultMaterial,
				 "TKR");	
  physiTKR = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+ACTSizeZ+
					     ACDTKRDistance+TKRSizeZ/2),
			       "TKR",		
			       logicTKR,
			       physiPayload,
			       false,	
			       0);	

  //                               
  // Anticoincidence Lateral Side (ACL)
  //
  
  solidACL1 = new G4Box("ACL1",			
			ACL1SizeX/2,ACL1SizeY/2,ACL1SizeZ/2); 
  
  logicACL1 = new G4LogicalVolume(solidACL1,ACDMaterial,"ACL");	
  
  physiACL1 = new G4PVPlacement(0, 
				G4ThreeVector(-PayloadSizeXY/2+ACL1SizeX/2,
					      -PayloadSizeXY/2+ACL1SizeY/2,
					      -PayloadSizeZ/2+
					      ACTSizeZ+ACL1SizeZ/2),
				"ACL1",		
				logicACL1,
				physiPayload,
				false,	
				0);	
  
  physiACL1 = new G4PVPlacement(0,
				G4ThreeVector(PayloadSizeXY/2-ACL1SizeX/2,
					      PayloadSizeXY/2-ACL1SizeY/2,
					      -PayloadSizeZ/2+
					      ACTSizeZ+ACL1SizeZ/2),
				"ACL1",		
				logicACL1,
				physiPayload,
				false,	
				1);	
  
  solidACL2 = new G4Box("ACL2",			
			ACL2SizeX/2,ACL2SizeY/2,ACL2SizeZ/2); 
  
  logicACL2 = new G4LogicalVolume(solidACL2,
				  ACDMaterial,
				  "ACL2");	
  
  
  physiACL2 = new G4PVPlacement(0, 
				G4ThreeVector(-PayloadSizeXY/2+ACL2SizeX/2,
					      PayloadSizeXY/2-ACL2SizeY/2,
					      -PayloadSizeZ/2+
					      ACTSizeZ+ACL2SizeZ/2),
				"ACL2",		
				logicACL2,
				physiPayload,
				false,	
				0);	
  
  physiACL2 = new G4PVPlacement(0, 
				G4ThreeVector(PayloadSizeXY/2-ACL2SizeX/2,
					      -PayloadSizeXY/2+ACL2SizeY/2,
					      -PayloadSizeZ/2+
					      ACTSizeZ+ACL2SizeZ/2),
				"ACL2",		
				logicACL2,
				physiPayload,
				false,	
				1);	
  
      
  // Tracker Structure (Plane + Converter + TKRDetectorX + TKRDetectorY)
  
  solidPlane = new G4Box("Plane",			
			 TKRSizeXY/2,TKRSizeXY/2,TKRSupportThickness/2); 
  
  logicPlane = new G4LogicalVolume(solidPlane,
				   defaultMaterial, 
				   "Plane");	
  

  solidTKRDetector = new G4Box
    ("TKRDetector",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  logicTKRDetector = new G4LogicalVolume(solidTKRDetector,
					 TKRMaterial, 
					 "TKRDetector");	
  
  
  solidConverter = new G4Box
    ("Converter",TKRSizeXY/2,TKRSizeXY/2,ConverterThickness/2); 
  
  logicConverter = new G4LogicalVolume(solidConverter,
				       ConverterMaterial, 
				       "Converter");	
  
  G4int i=0;
  
  for (i = 0; i < NbOfTKRLayers; i++)
    {

      physiPlane =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  (i)*TKRLayerDistance +
					  TKRSupportThickness/2),
			  "Plane",		
			  logicPlane,
			  physiTKR,
			  false,	
			  i);	

      physiConverter = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  ConverterThickness/2+
					  TKRSupportThickness+
					  (i)*TKRLayerDistance),
			  "Converter",		
			  logicConverter,
			  physiTKR,
			  false,	
			  i);
      
      physiTKRDetector = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  TKRSiliconThickness/2 +
					  ConverterThickness  +
					  TKRSupportThickness +
					  (i)*TKRLayerDistance),
			  "TKRDetector",		
			  logicTKRDetector,
			  physiTKR,
			  false,	
			  i);
      

      
      
      
    }
 


  G4VSolid * solidTKRActiveTile = new
    G4Box("Active Tile",TKRActiveTileXY/2,TKRActiveTileXY/2,TKRActiveTileZ/2);
  
  G4LogicalVolume* logicTKRActiveTile = 
    new G4LogicalVolume(solidTKRActiveTile, TKRMaterial,"Active Tile",0,0,0);
  
  
  G4int j=0;
  G4int k=0;
  
  G4VPhysicalVolume* physiTKRActiveTile = 0;
  
  G4double x=0.;
  G4double y=0.;
  G4double z=0.;
  
  for (i=0;i< NbOfTKRTiles; i++)
    { 
      for (j=0;j< NbOfTKRTiles; j++)
	{
	  k = i*NbOfTKRTiles + j;
	      
	  
	  x = -TKRSizeXY/2+TKRTileDistance+SiliconGuardRing+
	    TKRActiveTileXY/2+(i)*((2*SiliconGuardRing)+
				   TKRTileDistance+TKRActiveTileXY);
	  y = -TKRSizeXY/2+TKRTileDistance+SiliconGuardRing+
	    TKRActiveTileXY/2+(j)*((2*SiliconGuardRing)+TKRTileDistance+
				   TKRActiveTileXY);
	  z = 0.;
	  
	  	

	      
	  physiTKRActiveTile =
	    new G4PVPlacement(0,G4ThreeVector(x,y,z),
			      "Active Tile",logicTKRActiveTile,
			      physiTKRDetector,false, k);	
	  
	}
    }


      
  //                               
  // Sensitive Detectors: TKRActiveTile
  //
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(!trackerSD)
    {
      trackerSD = new Tst20TrackerSD("TrackerSD",this);
      SDman->AddNewDetector( trackerSD );		
    }
  G4String ROgeometryName = "TrackerROGeom";
  G4VReadOutGeometry* trackerRO = 
    new Tst20TrackerROGeometry(ROgeometryName, this);
  
  trackerRO->BuildROGeometry();

  trackerSD->SetROgeometry(trackerRO);
  if (logicTKRActiveTile)
    logicTKRActiveTile->SetSensitiveDetector(trackerSD); // sensitive tile

  
  if(!anticoincidenceSD)
    {
      anticoincidenceSD = new 
	Tst20AnticoincidenceSD("AnticoincidenceSD",this);
      SDman->AddNewDetector( anticoincidenceSD );		
    }
  
    
  if (logicACT)
    logicACT->SetSensitiveDetector(anticoincidenceSD); // top ACD
  if (logicACB)
    logicACB->SetSensitiveDetector(anticoincidenceSD); // botton ACD               
  if (logicACL1)
    logicACL1->SetSensitiveDetector(anticoincidenceSD); // lateral ACD
  if (logicACL2)
    logicACL2->SetSensitiveDetector(anticoincidenceSD); // lateral ACD
  
  //
  //always return the physical World
  //

  PrintPayloadParameters();
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::PrintPayloadParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The Tracker is " << NbOfTKRLayers << " layers of:  "
         << ConverterThickness/mm << "mm of " << ConverterMaterial->GetName() 
         << "\n------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::SetConverterMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
    {
      ConverterMaterial = pttoMaterial;
      logicConverter->SetMaterial(pttoMaterial); 
      PrintPayloadParameters();
    }             
}

void Tst20DetectorConstruction::SetConverterThickness(G4double val)
{
  ConverterThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::SetTKRSiliconThickness(G4double val)
{
  TKRSiliconThickness = val;
}  


void Tst20DetectorConstruction::SetTKRSiliconPitch(G4double val)
{
  TKRSiliconPitch = val;
}  


void Tst20DetectorConstruction::SetTKRTileSizeXY(G4double val)
{
  TKRSiliconTileXY = val;
}  


void Tst20DetectorConstruction::SetNbOfTKRLayers(G4int val)
{
  NbOfTKRLayers = val;
}


void Tst20DetectorConstruction::SetNbOfTKRTiles(G4int val)
{
  NbOfTKRTiles = val;
}

void Tst20DetectorConstruction::SetTKRLayerDistance(G4double val)
{
  TKRLayerDistance = val;
}

void Tst20DetectorConstruction::SetTKRTileDistance(G4double val)
{
  TKRTileDistance = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::SetACDThickness(G4double val)
{
  ACDThickness = val;
}

void Tst20DetectorConstruction::SetACDMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
    {
      ACDMaterial = pttoMaterial;
      logicACT->SetMaterial(pttoMaterial);
      logicACB->SetMaterial(pttoMaterial); 
      logicACL1->SetMaterial(pttoMaterial); 
      logicACL2->SetMaterial(pttoMaterial); 
      PrintPayloadParameters();
    }             
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if(magField) delete magField;		//delete the existing magn field
  
  if(fieldValue!=0.)			// create a new one if non nul
    { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    } else {
      magField = 0;
      fieldMgr->SetDetectorField(magField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void Tst20DetectorConstruction::UpdateGeometry()
{
  //  delete payloadSD;
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructPayload());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

















