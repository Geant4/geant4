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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelDetectorConstruction  ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelDetectorMessenger.hh"

#include "GammaRayTelTrackerSD.hh"

#include "GammaRayTelAnticoincidenceSD.hh"
#include "GammaRayTelCalorimeterSD.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UImanager.hh"


#include "G4RegionStore.hh"

G4ThreadLocal G4GlobalMagFieldMessenger* 
 GammaRayTelDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorConstruction::GammaRayTelDetectorConstruction()
  :solidWorld(0),logicWorld(0),physiWorld(0),
   solidPayload(0),logicPayload(0),physiPayload(0),
   solidTKR(0),logicTKR(0),physiTKR(0),
   solidCAL(0),logicCAL(0),physiCAL(0),
   solidACT(0),logicACT(0),physiACT(0),
   solidACL1(0),logicACL1(0),physiACL1(0),
   solidACL2(0),logicACL2(0),physiACL2(0),
   solidTKRDetectorX(0),logicTKRDetectorX(0),physiTKRDetectorX(0),
   solidTKRDetectorY(0),logicTKRDetectorY(0),physiTKRDetectorY(0),
   solidCALLayerX(0),logicCALLayerX(0),physiCALLayerX(0),
   solidCALLayerY(0),logicCALLayerY(0),physiCALLayerY(0),
   solidCALDetectorX(0),logicCALDetectorX(0),physiCALDetectorX(0),
   solidCALDetectorY(0),logicCALDetectorY(0),physiCALDetectorY(0),
   solidPlane(0),logicPlane(0),physiPlane(0),
   solidConverter(0),logicConverter(0),physiConverter(0),
   logicTKRStripX(0),logicTKRStripY(0)
   // aTKRRegion(0), aCALRegion(0)
{
  // default parameter values of the payload
  
  ConverterThickness = 300.*micrometer;
  TKRSiliconThickness = 400.*micrometer;
  TKRSiliconTileXY = 9.*cm;
  TKRSiliconPitch = 200.*micrometer; 
  TKRLayerDistance = 3.*cm;
  SiliconGuardRing = 1.5*mm;
  TKRViewsDistance = 1.*mm;
  NbOfTKRLayers = 15;
  NbOfTKRTiles = 4;
  CALBarThickness = 1.5*cm;
  NbOfCALBars = 12;
  NbOfCALLayers = 5;
  ACDThickness = 1.*cm;
  NbOfACDTopTiles = 1;
  NbOfACDLateralTiles = 2;

  TilesSeparation = 100.*micrometer;
  ACDTKRDistance = 5.*cm;
  CALTKRDistance = 1.5*cm;

  //Initialize thread-local sensitive detectors
  trackerSD.Put(0);
  calorimeterSD.Put(0);
  anticoincidenceSD.Put(0);

  ComputePayloadParameters();

  // create commands for interactive definition of the payload
  detectorMessenger = new GammaRayTelDetectorMessenger(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelDetectorConstruction::~GammaRayTelDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* GammaRayTelDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructPayload();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::DefineMaterials()
{ 
  
  G4String name, symbol;    
  G4double a, z, density;            
  
  G4int ncomponents, natoms;
  G4double fractionmass;
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

  a = 26.98*g/mole;
  G4Element* Alumin = new G4Element(name="Aluminum"  ,symbol="Al" , z= 13., a);

  a = 28.09*g/mole;
  G4Element* Silicon = new G4Element(name="Silicon", symbol="Si", z=14., a);
  
  a= 55.845*g/mole;
  G4Element* Iron = new G4Element(name="Iron", symbol="Fe", z=26.,a);

  a = 126.904*g/mole;
  G4Element* I  = new G4Element(name="Iodine"  ,symbol="I" , z= 53., a);
  
  a = 132.905*g/mole;
  G4Element* Cs  = new G4Element(name="Cesium"  ,symbol="Cs" , z= 55., a);

  a = 207.19*g/mole;
  G4Element* Lead = new G4Element(name="Lead", symbol="Pb", z=82., a);

  //
  // define simple materials
  //

  density = 19.3*g/cm3;
  a = 183.84*g/mole;
  G4Material* W = new G4Material(name="Tungsten", z=74., a, density);
  
      
  //
  // define a material from elements.   case 1: chemical molecule
  //
  
  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  
  density = 4.53*g/cm3;
  G4Material* CsI = new G4Material(name="CesiumIodide", density, ncomponents=2);
  CsI->AddElement(Cs, natoms=5);
  CsI->AddElement(I, natoms=5);
  
  //
  // define a material from elements.   case 2: mixture by fractional mass
  //
  
  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);

  density = 2.700*g/cm3;
  G4Material* Al = new G4Material(name="Aluminum", density, ncomponents=1);
  Al->AddElement(Alumin, fractionmass=1.);

  density = 2.333*g/cm3;  
  G4Material* Si = new G4Material(name="Silicon", density, ncomponents=1);
  Si->AddElement(Silicon, fractionmass=1.);
  
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Iron", density, ncomponents=1);
  Fe->AddElement(Iron, fractionmass=1.);
  
  density = 11.35*g/cm3;
  G4Material* Pb = new G4Material(name="Lead", density, ncomponents=1);
  Pb->AddElement(Lead, fractionmass=1.);

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

  ConverterMaterial = W;
  defaultMaterial  = vacuum;
  ACDMaterial = Sci;
  CALMaterial = CsI;
  TKRMaterial = Si;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* GammaRayTelDetectorConstruction::ConstructPayload()
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
  
  /* solidPayload=0; logicPayload=0; physiPayload=0;
  solidTKR=0;logicTKR=0;physiTKR=0;
  solidCAL=0;logicCAL=0;physiCAL=0;
  solidACT=0;logicACT=0;physiACT=0;
  solidACL1=0;logicACL1=0;physiACL1=0; 
  solidACL2=0;logicACL2=0;physiACL2=0; 
  solidConverter=0;logicConverter=0;physiConverter=0; 
  solidTKRDetectorX=0;logicTKRDetectorX=0;
  solidTKRDetectorY=0;logicTKRDetectorY=0;
  physiTKRDetectorX=0;physiTKRDetectorY=0;
  solidCALDetectorX=0;logicCALDetectorX=0;physiCALDetectorX=0;
  solidCALDetectorY=0;logicCALDetectorY=0;physiCALDetectorY=0;
  solidPlane=0;logicPlane=0;physiPlane=0;
  aCALRegion=0; aTKRRegion=0;
  */
  //
  // Payload
  //
  
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
  // Calorimeter (CAL)
  //
  
  solidCAL = new G4Box("CAL",			
		       CALSizeXY/2,CALSizeXY/2,CALSizeZ/2); 
  
  logicCAL = new G4LogicalVolume(solidCAL,
				 defaultMaterial,
				 "CAL");	
  physiCAL = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+CALSizeZ/2),
			       "CAL",		
			       logicCAL,
			       physiPayload,
			       false,	
			       0);	
  //                                 
  // Tracker (TKR)
      //
  
  solidTKR = new G4Box("TKR",			
			   TKRSizeXY/2,TKRSizeXY/2,TKRSizeZ/2); 
  
  logicTKR = new G4LogicalVolume(solidTKR,
				 defaultMaterial,
				 "TKR");	
  physiTKR = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+CALSizeZ+
					     CALTKRDistance+TKRSizeZ/2),
			       "TKR",		
			       logicTKR,
			       physiPayload,
			       false,	
			       0);	
  
  
  //                               
  // Anticoincidence Top (ACT)
  //
  
  solidACT = new G4Box("ACT",			
		       ACTSizeXY/2,ACTSizeXY/2,ACTSizeZ/2); 
  
  logicACT = new G4LogicalVolume(solidACT,ACDMaterial,"ACT");
  
  physiACT = new G4PVPlacement(0,		
			       G4ThreeVector(0,0,
					     -PayloadSizeZ/2+CALSizeZ+
					     CALTKRDistance+TKRSizeZ+
					     ACDTKRDistance+ACTSizeZ/2),
			       "ACT",		
			       logicACT,
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
					      -PayloadSizeZ/2+ACL1SizeZ/2),
				"ACL1",		
				logicACL1,
				physiPayload,
				false,	
				0);	
  
  physiACL1 = new G4PVPlacement(0,
				G4ThreeVector(PayloadSizeXY/2-ACL1SizeX/2,
					      PayloadSizeXY/2-ACL1SizeY/2,
					      -PayloadSizeZ/2+ACL1SizeZ/2),
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
					      -PayloadSizeZ/2+ACL2SizeZ/2),
				"ACL2",		
				logicACL2,
				physiPayload,
				false,	
				0);	
  
  physiACL2 = new G4PVPlacement(0, 
				G4ThreeVector(PayloadSizeXY/2-ACL2SizeX/2,
					      -PayloadSizeXY/2+ACL2SizeY/2,
					      -PayloadSizeZ/2+ACL2SizeZ/2),
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
  
  solidTKRDetectorY = new G4Box
    ("TKRDetectorY",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  logicTKRDetectorY = new G4LogicalVolume(solidTKRDetectorY,
					  TKRMaterial, 
					  "TKRDetector Y");	
  
  
  solidTKRDetectorX = new G4Box
    ("TKRDetectorX",TKRSizeXY/2,TKRSizeXY/2,TKRSiliconThickness/2); 
  
  logicTKRDetectorX = new G4LogicalVolume(solidTKRDetectorX,
					  TKRMaterial, 
					  "TKRDetector X");	
  
  
  solidConverter = new G4Box
    ("Converter",TKRSizeXY/2,TKRSizeXY/2,ConverterThickness/2); 
  
  logicConverter = new G4LogicalVolume(solidConverter,
				       ConverterMaterial, 
				       "Converter");	
  
  G4int i=0;
  
  for (i = 0; i < NbOfTKRLayers; i++)
    {
      
      physiTKRDetectorY = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,-TKRSizeZ/2
					  +TKRSiliconThickness/2 
					  +(i)*TKRLayerDistance),
			  "TKRDetectorY",		
			  logicTKRDetectorY,
			  physiTKR,
			  false,	
			  i);
      
      physiTKRDetectorX = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  TKRSiliconThickness/2 +
					  TKRViewsDistance+
					  TKRSiliconThickness+
					  (i)*TKRLayerDistance),
			  "TKRDetectorX",		
			  logicTKRDetectorX,
			  physiTKR,
			  false,	
			  i);
      
      
      physiConverter = 
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  2*TKRSiliconThickness +
					  TKRViewsDistance+
					  ConverterThickness/2+
					  (i)*TKRLayerDistance),
			  "Converter",		
			  logicConverter,
			  physiTKR,
			  false,	
			  i);

            
      
      physiPlane =
	new G4PVPlacement(0,G4ThreeVector(0.,0.,
					  -TKRSizeZ/2+
					  2*TKRSiliconThickness +
					  TKRViewsDistance+
					  ConverterThickness+
					  TKRSupportThickness/2+
					  (i)*TKRLayerDistance),
			  "Plane",		
			  logicPlane,
			  physiTKR,
			  false,	
			  i);	
      
    }
  
  
  
  G4VSolid * solidTKRActiveTileX = new
    G4Box("Active Tile X", TKRActiveTileXY/2,
	  TKRActiveTileXY/2,TKRActiveTileZ/2);
  
  
  G4VSolid * solidTKRActiveTileY = new
    G4Box("Active Tile Y", TKRActiveTileXY/2,
	  TKRActiveTileXY/2,TKRActiveTileZ/2);
  
  
  G4LogicalVolume* logicTKRActiveTileX = 
    new G4LogicalVolume(solidTKRActiveTileX, TKRMaterial,
			"Active Tile X",0,0,0);
  
  
  G4LogicalVolume* logicTKRActiveTileY = 
    new G4LogicalVolume(solidTKRActiveTileY, TKRMaterial,
			"Active Tile Y",0,0,0);
  
  


  G4int j=0;
  G4int k=0;
  
  G4double x=0.;
  G4double y=0.;
  G4double z=0.;
  
  for (i=0;i< NbOfTKRTiles; i++)
    { 
      for (j=0;j< NbOfTKRTiles; j++)
	{
	  k = i*NbOfTKRTiles + j;
	  
	  
	  x = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+
	    TKRActiveTileXY/2+(i)*((2*SiliconGuardRing)+
				   TilesSeparation+TKRActiveTileXY);
	  y = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+
	    TKRActiveTileXY/2+(j)*((2*SiliconGuardRing)+TilesSeparation+
				   TKRActiveTileXY);
	  z = 0.;
	  
	  new G4PVPlacement(0,
			    G4ThreeVector(x,y,z),
			    logicTKRActiveTileY,
			    "Active Tile Y",		
			    logicTKRDetectorY,
			    false,	
			    k);
	  
	  
	  x = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+
	    TKRActiveTileXY/2+(j)*((2*SiliconGuardRing)+
				   TilesSeparation+TKRActiveTileXY);
	  y = -TKRSizeXY/2+TilesSeparation+SiliconGuardRing+
	    TKRActiveTileXY/2+(i)*((2*SiliconGuardRing)+
				   TilesSeparation+TKRActiveTileXY);
	  z = 0.;
	      
	  new G4PVPlacement(0,
			    G4ThreeVector(x,y,z),
			    logicTKRActiveTileX,
			    "Active Tile X",
			    logicTKRDetectorX,
			    false,	
			    k);	
	  
	}
    }


  // STRIPS (not any more in the Readout Geometry)

  // Silicon Strips 
    
  /*
    G4double TKRXStripX=0.;
    G4double TKRYStripY=0.;
    G4double TKRYStripX=0.; 
    G4double TKRXStripY=0.;
  */

  TKRXStripX = TKRYStripY = TKRSiliconPitch;
  TKRYStripX = TKRXStripY = TKRActiveTileXY;
  TKRZStrip  = TKRSiliconThickness;
  
  
  G4VSolid* solidTKRStripX = new G4Box("Strip X",			
				       TKRXStripX/2,TKRYStripX/2,
				       TKRZStrip/2); 
  
  logicTKRStripX = 
    new G4LogicalVolume(solidTKRStripX,TKRMaterial,"Strip X",0,0,0);	 
  
		
  G4VSolid* solidTKRStripY = new G4Box("Strip Y",			
				       TKRXStripY/2,TKRYStripY/2,
				       TKRZStrip/2); 
  

  logicTKRStripY = 
    new G4LogicalVolume(solidTKRStripY,TKRMaterial,"Strip Y",0,0,0);	 
	

  for (i=0;i< NbOfTKRStrips; i++)
    {  
      new G4PVPlacement(0,
			G4ThreeVector(-TKRActiveTileXY/2 +TKRSiliconPitch/2 +
				      (i)*TKRSiliconPitch, 0., 0.),
			logicTKRStripX,
			"Strip X",		
			logicTKRActiveTileX,
			false,	
			i);	
      
      
      new G4PVPlacement(0,
			G4ThreeVector(0.,-TKRActiveTileXY/2 
				      +TKRSiliconPitch/2 +
				      (i)*TKRSiliconPitch, 0.),
			logicTKRStripY,
			"Strip Y",		
			logicTKRActiveTileY,
			false,	
			i);	
      
      
      
      

    }
  


  // Calorimeter Structure (CALLayerX + CALLayerY)
  
  
  solidCALLayerX = new G4Box("CALLayerX",			
			     CALSizeXY/2,CALSizeXY/2,CALBarThickness/2); 
  
  logicCALLayerX = new G4LogicalVolume(solidCALLayerX,
				       CALMaterial, 
				       "CALLayerX");	
  
  solidCALLayerY = new G4Box("CALLayerY",			
			     CALSizeXY/2,CALSizeXY/2,CALBarThickness/2); 
  
  logicCALLayerY = new G4LogicalVolume(solidCALLayerY,
				       CALMaterial, 
				       "CALLayerY");	
  
  for (i = 0; i < NbOfCALLayers; i++)
    {
      
      physiCALLayerY = 
	new G4PVPlacement(0,G4ThreeVector(0,0,
					  -CALSizeZ/2+
					  CALBarThickness/2 +
					  (i)*2*CALBarThickness),
			  "CALLayerY",		
			  logicCALLayerY,
			  physiCAL,
			  false,	
			  i);	
      
      physiCALLayerX = 
	new G4PVPlacement(0,G4ThreeVector(0,0,
					  -CALSizeZ/2+
					  CALBarThickness/2 + 
					  CALBarThickness +
					  (i)*2*CALBarThickness),
			  "CALLayerX",		
			  logicCALLayerX,
			  physiCAL,
			  false,	
			  i);	
      
    }
  
  // Calorimeter Structure (CALDetectorX + CALDetectorY)
  
  solidCALDetectorX = new G4Box("CALDetectorX",			
				CALBarX/2,CALBarY/2,CALBarThickness/2); 
  
  logicCALDetectorX = new G4LogicalVolume(solidCALDetectorX,
					  CALMaterial, 
					  "CALDetectorX");	
  
  solidCALDetectorY = new G4Box("CALDetectorY",			
				CALBarY/2,CALBarX/2,CALBarThickness/2); 
  
  logicCALDetectorY = new G4LogicalVolume(solidCALDetectorY,
					  CALMaterial, 
					  "CALDetectorY");	
  
  for (i = 0; i < NbOfCALBars; i++)
    {
      
      physiCALDetectorY = 
	new G4PVPlacement(0,
			  G4ThreeVector(-CALSizeXY/2+ CALBarY/2 +
					(i)*CALBarY, 0, 0),
			  logicCALDetectorY,
			  "CALDetectorY",		
			  logicCALLayerY,
			  false,	
			  i);	
      
      physiCALDetectorX = 
	new G4PVPlacement(0,
			  G4ThreeVector(0,-CALSizeXY/2+ CALBarY/2 +
					(i)*CALBarY, 0),
			  logicCALDetectorX,
			  "CALDetectorX",		
			  logicCALLayerX,
			  false,	
			  i);	
      
    }
  
  
  // Cuts by Regions 

  /*
    G4String regName[] = {"Calorimeter","Tracker"};
    if (aCALRegion) delete aCALRegion;
    
    aCALRegion = new G4Region(regName[0]);
    logicCAL->SetRegion(aCALRegion);
    aCALRegion->AddRootLogicalVolume(logicCAL);
    
    if (aTKRRegion) delete aTKRRegion;
    
    aTKRRegion = new G4Region(regName[1]);
    logicTKR->SetRegion(aTKRRegion);
    aTKRRegion->AddRootLogicalVolume(logicTKR);
  */


  //                                        
  // Visualization attributes
  //
  
  // Invisible Volume
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicPayload->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicTKR->SetVisAttributes(G4VisAttributes::GetInvisible());  
  logicTKRActiveTileX->SetVisAttributes(G4VisAttributes::GetInvisible());  
  logicTKRActiveTileY->SetVisAttributes(G4VisAttributes::GetInvisible());  
  logicPlane->SetVisAttributes(G4VisAttributes::GetInvisible());  
  logicConverter->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicCAL->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicCALLayerX->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicCALLayerY->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicTKRStripX->SetVisAttributes(G4VisAttributes::GetInvisible());
  logicTKRStripY->SetVisAttributes(G4VisAttributes::GetInvisible());  

  // Some visualization styles

  G4VisAttributes* VisAtt1= new G4VisAttributes(G4Colour(0.3,0.8,0.1));
  VisAtt1->SetVisibility(true);
  VisAtt1->SetForceSolid(TRUE);

  G4VisAttributes* VisAtt2= new G4VisAttributes(G4Colour(0.2,0.3,0.8));
  VisAtt2->SetVisibility(true);
  VisAtt2->SetForceSolid(FALSE);

  G4VisAttributes* VisAtt3= new G4VisAttributes(G4Colour(0.8,0.2,0.3));
  VisAtt3->SetVisibility(true);
  VisAtt3->SetForceWireframe(TRUE);
  
  // Visible Volumes

  logicCALDetectorX->SetVisAttributes(VisAtt1);
  logicCALDetectorY->SetVisAttributes(VisAtt1);
  logicTKRDetectorX->SetVisAttributes(VisAtt2);
  logicTKRDetectorY->SetVisAttributes(VisAtt2);
  logicACT->SetVisAttributes(VisAtt3);  
  logicACL1->SetVisAttributes(VisAtt3);  
  logicACL2->SetVisAttributes(VisAtt3);


  //
  //always return the physical World
  //
  PrintPayloadParameters();
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::ConstructSDandField()
{ 
  
  //
  // Sensitive Detectors - Tracker
  //                                
  if(trackerSD.Get()==0)
    {
      GammaRayTelTrackerSD* SD = new GammaRayTelTrackerSD("TrackerSD");
      trackerSD.Put(SD);
    }

  G4SDManager::GetSDMpointer()->AddNewDetector(trackerSD.Get());
  //Flags the strips as sensitive .
  if (logicTKRStripX)
    SetSensitiveDetector(logicTKRStripX,trackerSD.Get()); // ActiveStripX
  if (logicTKRStripY)
    SetSensitiveDetector(logicTKRStripY,trackerSD.Get()); // ActiveStripY
  

  //
  // Sensitive Detectors: Calorimeter
  // 
  if(calorimeterSD.Get()==0)
    {
      GammaRayTelCalorimeterSD* SD = new GammaRayTelCalorimeterSD("CalorimeterSD");
      calorimeterSD.Put(SD);
    }  
  G4SDManager::GetSDMpointer()->AddNewDetector(calorimeterSD.Get());
  if (logicCALDetectorX)
    SetSensitiveDetector(logicCALDetectorX,calorimeterSD.Get()); // BarX
  if (logicCALDetectorY)
    SetSensitiveDetector(logicCALDetectorY,calorimeterSD.Get()); // BarY

  //
  // Sensitive Detectors: Anticoincidence
  //

  if(anticoincidenceSD.Get()==0)
    {
      GammaRayTelAnticoincidenceSD* SD = new GammaRayTelAnticoincidenceSD
	("AnticoincidenceSD");
      anticoincidenceSD.Put(SD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(anticoincidenceSD.Get()); 
  if (logicACT)
    SetSensitiveDetector(logicACT,anticoincidenceSD.Get()); // ACD top
  if (logicACL1)
    SetSensitiveDetector(logicACL1,anticoincidenceSD.Get()); // ACD lateral side
  if (logicACL2)
    SetSensitiveDetector(logicACL2,anticoincidenceSD.Get()); // ACD lateral side

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

  return;
  
}

void GammaRayTelDetectorConstruction::PrintPayloadParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The Tracker is " << NbOfTKRLayers << " layers of:  "
         << ConverterThickness/mm << "mm of " << ConverterMaterial->GetName() 
         << "\n------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetConverterMaterial(G4String materialChoice)
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

void GammaRayTelDetectorConstruction::SetConverterThickness(G4double val)
{
  ConverterThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetTKRSiliconThickness(G4double val)
{
  TKRSiliconThickness = val;
}  


void GammaRayTelDetectorConstruction::SetTKRSiliconPitch(G4double val)
{
  TKRSiliconPitch = val;
}  


void GammaRayTelDetectorConstruction::SetTKRTileSizeXY(G4double val)
{
  TKRSiliconTileXY = val;
}  


void GammaRayTelDetectorConstruction::SetNbOfTKRLayers(G4int val)
{
  NbOfTKRLayers = val;
}


void GammaRayTelDetectorConstruction::SetNbOfTKRTiles(G4int val)
{
  NbOfTKRTiles = val;
}

void GammaRayTelDetectorConstruction::SetTKRLayerDistance(G4double val)
{
  TKRLayerDistance = val;
}

void GammaRayTelDetectorConstruction::SetTKRViewsDistance(G4double val)
{
  TKRViewsDistance = val;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetNbOfCALLayers(G4int val)
{
  NbOfCALLayers = val;
}

void GammaRayTelDetectorConstruction::SetNbOfCALBars(G4int val)
{
  NbOfCALBars = val;
}

void GammaRayTelDetectorConstruction::SetCALBarThickness(G4double val)
{
  CALBarThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetACDThickness(G4double val)
{
  ACDThickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::SetMagField(G4double fieldValue)
{
  // Just invoke manually the MT-safe command 
  // /globalField/setValue 
  // instantiated by the GlobalFieldMessenger
  std::stringstream sss;
  sss << "/globalField/setValue 0 0 " << fieldValue/tesla << " tesla";

  G4String command = sss.str();
  G4cout << "Going to execute: " << command << G4endl;

  G4UImanager* UImanager = G4UImanager::GetUIpointer();  
  UImanager->ApplyCommand(command);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void GammaRayTelDetectorConstruction::UpdateGeometry()
{
  //  delete payloadSD;
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructPayload());
  G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  G4RegionStore::GetInstance()->UpdateMaterialList(physiWorld);

  G4RunManager::GetRunManager()->ReinitializeGeometry();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

















