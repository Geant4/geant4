// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelDetectorConstruction.cc,v 1.1 2000-10-05 09:48:07 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelDetectorMessenger.hh"

#include "GammaRayTelPayloadSD.hh"

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

GammaRayTelDetectorConstruction::GammaRayTelDetectorConstruction()
  :solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
   solidPayload(NULL),logicPayload(NULL),physiPayload(NULL),
   solidTKR(NULL),logicTKR(NULL),physiTKR(NULL),
   solidCAL(NULL),logicCAL(NULL),physiCAL(NULL),
   solidACT(NULL),logicACT(NULL),physiACT(NULL),
   solidACL1(NULL),logicACL1(NULL),physiACL1(NULL),
   solidACL2(NULL),logicACL2(NULL),physiACL2(NULL),
   solidConverter(NULL),logicConverter(NULL),physiConverter(NULL),
   solidTKRDetector(NULL),logicTKRDetector(NULL),
   physiTKRDetectorX(NULL),physiTKRDetectorY(NULL),
   solidTKRActiveTile(NULL),logicTKRActiveTile(NULL),
   physiTKRActiveTileX(NULL), physiTKRActiveTileY(NULL),
   solidTKRStripX(NULL),logicTKRStripX(NULL),physiTKRStripX(NULL),
   solidTKRStripY(NULL),logicTKRStripY(NULL),physiTKRStripY(NULL),
   solidCALDetector(NULL),logicCALDetector(NULL),
   physiCALDetectorX(NULL),physiCALDetectorY(NULL),
   solidPlane(NULL),logicPlane(NULL),physiPlane(NULL)

{
  // default parameter values of the payload
  
  ConverterThickness = 300.*micrometer;
  TKRSiliconThickness = 400.*micrometer;
  TKRSiliconTileXY = 9.*cm;
  TKRSiliconPitch = 200.*micrometer; // come fare per le strip?
  TKRLayerDistance = 3.*cm;
  TKRViewsDistance = 1.*mm;
  NbOfTKRLayers = 15;
  NbOfTKRTiles = 4;
  CALBarThickness = 1.5*cm;
  NbOfCALBars = 12;
  NbOfCALLayers = 5;
  ACDThickness = 1.*cm;

  TilesSeparation = 100.*micrometer;
  ACDTKRDistance = 5.*cm;
  CALTKRDistance = 1.5*cm;

  ComputePayloadParameters();

  // create commands for interactive definition of the calorimeter  
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
  solidWorld = new G4Box("World",				     //its name
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);  //its size
  
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Payload
  //  

  solidPayload=NULL; logicPayload=NULL; physiPayload=NULL;
  solidTKR=NULL;logicTKR=NULL;physiTKR=NULL;
  solidCAL=NULL;logicCAL=NULL;physiCAL=NULL;
  solidACT=NULL;logicACT=NULL;physiACT=NULL;
  solidACL1=NULL;logicACL1=NULL;physiACL1=NULL; 
  solidACL2=NULL;logicACL2=NULL;physiACL2=NULL; 
  solidConverter=NULL;logicConverter=NULL;physiConverter=NULL; 
  solidTKRDetector=NULL;logicTKRDetector=NULL;
  physiTKRDetectorX=NULL;physiTKRDetectorY=NULL;
  solidTKRActiveTile=NULL;logicTKRActiveTile=NULL;
  physiTKRActiveTileX=NULL;physiTKRActiveTileY=NULL;
  solidTKRStripX=NULL;logicTKRStripX=NULL;physiTKRStripX=NULL;
  solidTKRStripY=NULL;logicTKRStripY=NULL;physiTKRStripY=NULL;

  solidCALDetector=NULL;logicCALDetector=NULL;
  physiCALDetectorX=NULL;physiCALDetectorY=NULL;
  solidPlane=NULL;logicPlane=NULL;physiPlane=NULL;
  
  if (PayloadSizeZ > 0.) 
    { 
      

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
                                                 ACDTKRDistance+ACTSizeZ),
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

      physiACL1 = new G4PVPlacement(0, // <-- da definire rotation matrix
				    G4ThreeVector(-PayloadSizeXY/2+ACL1SizeX/2,
						  -PayloadSizeXY/2+ACL1SizeY/2,
						  -PayloadSizeZ/2+ACL1SizeZ/2),
				    "ACL1",		
				    logicACL1,
				    physiPayload,
				    false,	
				    0);	

      physiACL1 = new G4PVPlacement(0, // <-- da definire rotation matrix
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
 

      physiACL2 = new G4PVPlacement(0, // <-- da definire rotation matrix
				    G4ThreeVector(-PayloadSizeXY/2+ACL2SizeX/2,
						  PayloadSizeXY/2-ACL2SizeY/2,
						  -PayloadSizeZ/2+ACL2SizeZ/2),
				    "ACL2",		
				    logicACL2,
				    physiPayload,
				    false,	
				    0);	
      
      physiACL2 = new G4PVPlacement(0, // <-- da definire rotation matrix
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
				       defaultMaterial, // (devo fare le fibre)
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
      G4double dist = 0;
      for (i = 0; i < NbOfTKRLayers; i++)
	{
          G4cout << i << " plane" << G4endl;
	  dist = -TKRSizeZ/2+TKRSiliconThickness/2 +(i)*TKRLayerDistance;
	  G4cout << dist/cm << " distance" << G4endl;
          G4cout << TKRSizeZ/cm << " tracker" << G4endl;
          G4cout << TKRLayerDistance/cm << " layer" << G4endl;


          physiTKRDetectorY = 
	    new G4PVPlacement(0,G4ThreeVector(0,0,-TKRSizeZ/2
					      +TKRSiliconThickness/2 
					      +(i)*TKRLayerDistance),
			      "TKRDetectorY",		
			      logicTKRDetector,
			      physiTKR,
			      false,	
			      i);	
	  
          physiTKRDetectorX = 
	    new G4PVPlacement(0,G4ThreeVector(0.,0.,
					      -TKRSizeZ/2+
					      TKRSiliconThickness/2 +
					      TKRViewsDistance+
					      (i)*TKRLayerDistance),
			      "TKRDetectorX",		
			      logicTKRDetector,
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
			      logicTKRDetector,
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
                                              TKRSupportThickness/2),
			      "Plane",		
			      logicPlane,
			      physiTKR,
			      false,	
			      i);	
	  

	}
      
      solidTKRActiveTile = new G4Box("Active Tile",			
				     TKRActiveTileXY/2,TKRActiveTileXY/2,
				     TKRActiveTileZ/2); 
      
      logicTKRActiveTile = new G4LogicalVolume(solidTKRActiveTile,
					       TKRMaterial, 
					       "Active Tile");	

      G4int j=0;
      G4int k=0;

      for (i=0;i< NbOfTKRTiles; i++)
	{ 
	  for (j=0;j< NbOfTKRTiles; j++)
            {
              k = i*NbOfTKRTiles + j;
	      physiTKRActiveTileX =
		new G4PVPlacement(0,
				  G4ThreeVector(-TKRSizeXY/2+TilesSeparation+
						SiliconGuardRing+
						TKRActiveTileXY/2+
						(j)*((2*SiliconGuardRing)+
						     TilesSeparation+
						     TKRActiveTileXY),
						-TKRSizeXY/2+TilesSeparation+
						SiliconGuardRing+
						TKRActiveTileXY/2+
						(i)*((2*SiliconGuardRing)+
						     TilesSeparation+
						     TKRActiveTileXY),
						0.),
			      "Active Tile X",		
			      logicTKRActiveTile,
			      physiTKRDetectorX,
			      false,	
			      k);	
    
	      physiTKRActiveTileY =
		new G4PVPlacement(0,
				  G4ThreeVector(-TKRSizeXY/2+TilesSeparation+
						SiliconGuardRing+
						TKRActiveTileXY/2+
						(j)*((2*SiliconGuardRing)+
						     TilesSeparation+
						     TKRActiveTileXY),
						-TKRSizeXY/2+TilesSeparation+
						SiliconGuardRing+
						TKRActiveTileXY/2+
						(i)*((2*SiliconGuardRing)+
						     TilesSeparation+
						     TKRActiveTileXY),
						0.),
			      "Active Tile Y",		
			      logicTKRActiveTile,
			      physiTKRDetectorY,
			      false,	
			      k);	
	    }
	}



      solidTKRStripX = new G4Box("Strip X",			
				     TKRXStripX/2,TKRYStripX/2,
				     TKRZStrip/2); 
      
      logicTKRStripX = new G4LogicalVolume(solidTKRStripX,
					       TKRMaterial, 
					       "Strip X");	

      physiTKRStripX = new G4PVReplica("Strip X",       
                                     logicTKRStripX,       
                                     physiTKRActiveTileX,  
                                     kXAxis,            
                                     NbOfTKRStrips,
                                     TKRXStripX); 

      solidTKRStripY = new G4Box("Strip Y",			
				     TKRXStripY/2,TKRYStripY/2,
				     TKRZStrip/2); 
      
      logicTKRStripY = new G4LogicalVolume(solidTKRStripY,
					       TKRMaterial, 
					       "Strip Y");	

      physiTKRStripY = new G4PVReplica("Strip Y",       
                                     logicTKRStripY,       
                                     physiTKRActiveTileY,  
                                     kYAxis,            
                                     NbOfTKRStrips,
                                     TKRYStripY); 

      // Calorimeter Structure (CALDetectorX + CALDetectorY)
      
   
      solidCALDetector = new G4Box("CALDetector",			
			     CALSizeXY/2,CALSizeXY/2,CALBarThickness/2); 
      
      logicCALDetector = new G4LogicalVolume(solidPlane,
					     CALMaterial, 
					     "CALDetector");	
      
      
      G4double dist1 = 0;
      G4double dist2 = 0;

      for (i = 0; i < NbOfCALLayers; i++)
	{
	  
	  
          dist1 = -CALSizeZ/2+CALBarThickness/2+(i)*2*CALBarThickness;
          dist2 = -CALSizeZ/2+CALBarThickness/2+(i)*2*CALBarThickness+                           CALBarThickness;
 
          G4cout << dist1/mm << "Distanza 1" << G4endl;
          G4cout << dist2/mm << "Distanza 2" << G4endl;
          G4cout << CALSizeZ/mm << "Distanza CAL" << G4endl;

	  physiCALDetectorY = 
	    new G4PVPlacement(0,G4ThreeVector(0,0,
					      -CALSizeZ/2+
					      CALBarThickness/2 +
					      (i)*2*CALBarThickness),
			      "CALDetectorY",		
			      logicCALDetector,
			      physiCAL,
			      false,	
			      i);	
	  
          physiCALDetectorX = 
	    new G4PVPlacement(0,G4ThreeVector(0,0,
					      -CALSizeZ/2+
                                              CALBarThickness/2 + 
                                              CALBarThickness +
					      (i)*2*CALBarThickness),
			      "CALDetectorX",		
			      logicCALDetector,
			      physiCAL,
			      false,	
			      i);	

	}


    }

  //                               
  // Sensitive Detectors: TKRDetector CALDetector ACT ACL1 ACL2
  //
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if(!payloadSD)
    {
      payloadSD = new GammaRayTelPayloadSD("PayloadSD",this);
      SDman->AddNewDetector( payloadSD );
    }

  //  if (logicTKRDetector)
  //  logicTKRDetector->SetSensitiveDetector(payloadSD);

  if (logicTKRStripX)
    logicTKRStripX->SetSensitiveDetector(payloadSD);
  if (logicTKRStripY)
    logicTKRStripY->SetSensitiveDetector(payloadSD);

  /*
  if (logicCALDetector)
    logicCALDetector->SetSensitiveDetector(payloadSD);
  
  if (logicACT)
    logicACT->SetSensitiveDetector(payloadSD);
  if (logicACL1)
    logicACL1->SetSensitiveDetector(payloadSD);
  if (logicACL2)
    logicACL2->SetSensitiveDetector(payloadSD);
  */
  
  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicPayload->SetVisAttributes(simpleBoxVisAtt);
  
  //
  //always return the physical World
  //
  PrintPayloadParameters();
  return physiWorld;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelDetectorConstruction::PrintPayloadParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The Payload is " << NbOfTKRLayers << " layers of:  "
         << ConverterThickness/mm << "mm of " << ConverterMaterial->GetName() 
         << "\n------------------------------------------------------------\n";
}

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
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr 
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if(magField) delete magField;		//delete the existing magn field
  
  if(fieldValue!=0.)			// create a new one if non nul
    { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
    } else {
      magField = NULL;
      fieldMgr->SetDetectorField(magField);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void GammaRayTelDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructPayload());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




