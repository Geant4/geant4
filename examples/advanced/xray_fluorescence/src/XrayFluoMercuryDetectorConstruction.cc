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
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// 08 Sep 2003 Alfonso Mantero Created
// -------------------------------------------------------------------

#include "XrayFluoMercuryDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorMessenger.hh"
#include "XrayFluoSD.hh"
#include "XrayFluoNistMaterials.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoMercuryDetectorConstruction::XrayFluoMercuryDetectorConstruction()
  : detectorType(0),mercuryGranularity(false), DeviceSizeX(0),
    DeviceSizeY(0),DeviceThickness(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidScreen(0),logicScreen(0),physiScreen(0),
    solidMercury (0),logicMercury(0),physiMercury (0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    screenMaterial(0),OhmicPosMaterial(0), OhmicNegMaterial(0),
    pixelMaterial(0),mercuryMaterial(0),
    defaultMaterial(0),HPGeSD(0)
  
{ 
  materials = XrayFluoNistMaterials::GetInstance();

  DefineDefaultMaterials();

  NbOfPixelRows     =  1; // should be 1
  NbOfPixelColumns  =  1; // should be 1
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
  PixelSizeXY       = std::sqrt(40.) * mm *0.5e6; // should be std::sqrt(40) * mm
  PixelThickness = 3.5 * mm * 1e6; //should be 3.5 mm
  
  G4cout << "PixelThickness(mm): "<< PixelThickness/mm << G4endl;
  G4cout << "PixelSizeXY(cm): "<< PixelSizeXY/cm << G4endl;
  
  ContactSizeXY  = std::sqrt(40.) * mm * 0.5e6; //should be the same as PixelSize or lower 

  mercuryDia = 2 * 4880 * km ;
  sunDia =  1390000 * km ;
  mercurySunDistance = 57910000 * km ;

  
  OhmicNegThickness = 0.005*mm *0.5e6 ;
  OhmicPosThickness = 0.005*mm *0.5e6 ;
  
  screenThickness = 5 * mm *0.5e6;
  
  ThetaHPGe = 135. * deg ;
  PhiHPGe = 225. * deg  ;
  
  
  distDe = (mercuryDia/2 + 400 * km);
  
  distScreen = distDe + (screenThickness+PixelThickness)/2+OhmicPosThickness ;

  distOptic = distDe - 1.*m * 1e5;//!!!
  
  opticThickness = 1.* cm *0.5e6;
  opticDia = 21. * cm *0.5e6;  
  opticAperture = 1. * deg;

  PixelCopyNb=0;
  grainCopyNb=0;
  G4String defaultDetectorType = "sili";
  ComputeApparateParameters();
  SetDetectorType(defaultDetectorType);
  
  // create commands for interactive definition of the apparate
  
  detectorMessenger = new XrayFluoMercuryDetectorMessenger(this);
  G4cout << "XrayFluoMercuryDetectorConstruction created" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoMercuryDetectorConstruction* XrayFluoMercuryDetectorConstruction::instance = 0;

XrayFluoMercuryDetectorConstruction* XrayFluoMercuryDetectorConstruction::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoMercuryDetectorConstruction;
      
    }
  return instance;
}

void XrayFluoMercuryDetectorConstruction::SetDetectorType(G4String type)
{
  
  if (type=="sili")
    {
      detectorType = XrayFluoSiLiDetectorType::GetInstance();
    }
  else if (type=="hpge")
    {
      detectorType = XrayFluoHPGeDetectorType::GetInstance();
    }
  else 
    {
      G4ExceptionDescription execp;
      execp <<  type + "detector type unknown";
      G4Exception("XrayFluoMercuryDetectorConstruction::SetDetectorType()","example-xray_fluorescence05",
	  FatalException, execp);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoVDetectorType* XrayFluoMercuryDetectorConstruction::GetDetectorType() const
{
  return detectorType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryDetectorConstruction::~XrayFluoMercuryDetectorConstruction()
  
{ 
  delete detectorMessenger;
  delete detectorType;
  G4cout << "XrayFluoMercuryDetectorConstruction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoMercuryDetectorConstruction::Construct()
{
  return ConstructApparate();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorConstruction::DefineDefaultMaterials()
{
  
  
  //define materials of the apparate
  
  mercuryMaterial = materials->GetMaterial("Anorthosite");
  screenMaterial = materials->GetMaterial("G4_Pb");
  pixelMaterial = materials->GetMaterial("G4_Si");
  OhmicPosMaterial = materials->GetMaterial("G4_Cu");
  OhmicNegMaterial = materials->GetMaterial("G4_Pb");
  defaultMaterial = materials->GetMaterial("G4_Galactic");
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* XrayFluoMercuryDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition 
  
  ComputeApparateParameters();
  
  //world
  
  solidWorld = new G4Box("World",	      		        //its name
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);	//its size
  
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
  physiWorld = new G4PVPlacement(0,			//no rotation
				 G4ThreeVector(),	//at (0,0,0)
				 "World",		//its name
				 logicWorld,		//its logical volume
				 0,			//its mother  volume
				 false,			//no boolean operation
				 0);			//copy number
  
  //detector
  
  solidHPGe = 0;  physiHPGe = 0;  logicHPGe=0;
  solidPixel=0; logicPixel=0; physiPixel=0;
  
  if (DeviceThickness > 0.)  
    {
      solidHPGe = new G4Box("HPGeDetector",		//its name
			    DeviceSizeX/2,DeviceSizeY/2,DeviceThickness/2);//size
      
      
      logicHPGe = new G4LogicalVolume(solidHPGe,	//its solid
				      defaultMaterial,	//its material 
				      "HPGeDetector");	//its name
      
      zRotPhiHPGe.rotateX(PhiHPGe);
      G4double x,y,z;
      
      z = distDe * std::cos(ThetaHPGe);
      y = distScreen * std::sin(ThetaHPGe);
      x = 0.*cm;
      
      physiHPGe = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)), 
				    "HPGeDetector",	//its name
				    logicHPGe,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
    }
  // Pixel
  
  
  
  
  for ( G4int j=0; j < NbOfPixelColumns ; j++ )
    { for ( G4int i=0; i < NbOfPixelRows ; i++ )
      { 
	solidPixel=0; logicPixel=0;   physiPixel=0;
	if (PixelThickness > 0.)
	  solidPixel = new G4Box("Pixel",			
				 PixelSizeXY/2,PixelSizeXY/2, PixelThickness/2);
	
	logicPixel = new G4LogicalVolume(solidPixel,	
					 pixelMaterial,	//its material
					 "Pixel");	        //its name
	
	/*
	  zRotPhiHPGe.rotateX(PhiHPGe);
	  G4double x,y,z;
	  z = distDe * std::cos(ThetaHPGe);
	  y =distDe * std::sin(ThetaHPGe);
	  x = 0.*cm;*/ 
	physiPixel = new G4PVPlacement(0,	       
				       G4ThreeVector(0,
						     i*PixelSizeXY, 
						     j*PixelSizeXY ),
				       "Pixel",  
				       logicPixel,	 //its logical volume
				       physiHPGe, //its mother  volume
				       false,	 //no boolean operation
				       PixelCopyNb);//copy number
	
	
	
	
	
	
	// OhmicNeg
	
	solidOhmicNeg=0; logicOhmicNeg=0; physiOhmicNeg=0;  
	
	if (OhmicNegThickness > 0.) 
	  { solidOhmicNeg = new G4Box("OhmicNeg",		//its name
				      PixelSizeXY/2,PixelSizeXY/2,OhmicNegThickness/2); 
	  
	  logicOhmicNeg = new G4LogicalVolume(solidOhmicNeg,    //its solid
					      OhmicNegMaterial, //its material
					      "OhmicNeg");      //its name
	  
	  physiOhmicNeg = new G4PVPlacement(0,
					    G4ThreeVector
					    (0.,
					     0.,
					     (PixelThickness+OhmicNegThickness)/2),
					    "OhmicNeg",        //its name
					    logicOhmicNeg,     //its logical volume
					    physiHPGe,        //its mother
					    false,             //no boulean operat
					    PixelCopyNb);                //copy number
	  
	  }
	// OhmicPos
	solidOhmicPos=0; logicOhmicPos=0; physiOhmicPos=0;  
	
	if (OhmicPosThickness > 0.) 
	  { solidOhmicPos = new G4Box("OhmicPos",		//its name
				      PixelSizeXY/2,PixelSizeXY/2,OhmicPosThickness/2); 
	  
	  logicOhmicPos = new G4LogicalVolume(solidOhmicPos,    //its solid
					      OhmicPosMaterial, //its material
					      "OhmicPos");      //its name
	  
	  physiOhmicPos = new G4PVPlacement(0,	
					    G4ThreeVector(0.,
							  0.,
							  (-PixelThickness-OhmicPosThickness)/2),  
					    "OhmicPos",  
					    logicOhmicPos,
					    physiHPGe,  
					    false,     
					    PixelCopyNb); 
	  
	  }
	
	PixelCopyNb += PixelCopyNb; 
	G4cout << "PixelCopyNb: " << PixelCopyNb << G4endl;
      }
    
    }
  
  // Optics

  if (DeviceThickness > 0.)  
    {
      solidOptic = new G4Tubs("DetectorOptic",		//its name
			      0.,opticDia/2, opticThickness, 0.,2.*pi);//size
      
      
      logicOptic = new G4LogicalVolume(solidOptic,	//its solid
				       defaultMaterial,	//its material 
				       "DetectorOptic");	//its name
      
      //zRotPhiHPGe.rotateX(PhiHPGe);
      G4double x,y,z;
      z = distOptic * std::cos(ThetaHPGe);
      y = distOptic * std::sin(ThetaHPGe);
      x = 0.*cm;
      physiOptic = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)), 
				      "DetectorOptic",	//its name
				      logicOptic,	//its logical volume
				      physiWorld,	//its mother  volume
				      false,		//no boolean operation
				      0);		//copy number
    }
  

  // Screen
  
  if (DeviceThickness > 0.)  
    {
      solidScreen = new G4Box("DetectorScreen",		//its name
			      screenSizeXY/2,screenSizeXY/2,screenThickness/2);//size
      
      
      logicScreen = new G4LogicalVolume(solidScreen,	//its solid
					defaultMaterial,	//its material 
					"DetectorScreen");	//its name
      
      //zRotPhiHPGe.rotateX(PhiHPGe);
      G4double x,y,z;
      G4cout << "distScreen: "<< distScreen/m <<G4endl;
      z = distScreen * std::cos(ThetaHPGe);
      y = distScreen * std::sin(ThetaHPGe);
      x = 0.*cm;
      physiScreen = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)), 
				      "DetectorScreen",	//its name
				      logicScreen,	//its logical volume
				      physiWorld,	//its mother  volume
				      false,		//no boolean operation
				      0);		//copy number
    }
  
  //Mercury
  
  
  solidMercury=0;  logicMercury=0;  physiMercury=0;
  if (mercuryDia > 0.)  
    {


        



      solidMercury = new G4Sphere("Mercury",0.,mercuryDia/2., 0., twopi, 0., pi);
            
      logicMercury= new G4LogicalVolume(solidMercury,	//its solid
					mercuryMaterial,	//its material
					"Mercury");	//its name
      
      physiMercury = new G4PVPlacement(0,			//no rotation
				       G4ThreeVector(),	//at (0,0,0)
				       "Mercury",	//its name
				       logicMercury,	//its logical volume
				       physiWorld,	//its mother  volume
				       false,		//no boolean operation
				       0);		//copy number
      
    }  
  
   
  // Visualization attributes
  

  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255. ,255/255. ,51/255. ));
  G4VisAttributes * red= new G4VisAttributes( G4Colour(255/255. , 0/255. , 0/255. ));
  G4VisAttributes * blue= new G4VisAttributes( G4Colour(0/255. , 0/255. ,  255/255. ));
  G4VisAttributes * grayc= new G4VisAttributes( G4Colour(128/255. , 128/255. ,  128/255. ));
  G4VisAttributes * darkGray= new G4VisAttributes( G4Colour(95/255. , 95/255. ,  95/255. ));
  //G4VisAttributes * green= new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  red->SetVisibility(true);
  red->SetForceSolid(true);
  blue->SetVisibility(true);
  grayc->SetVisibility(true);
  grayc->SetForceSolid(true);
  simpleBoxVisAtt->SetVisibility(true);

  //logicWorld->SetVisAttributes (simpleBoxVisAtt);
  
  logicPixel->SetVisAttributes(red);
  logicHPGe->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  logicMercury->SetVisAttributes(darkGray);
  

  logicScreen->SetVisAttributes(red);
  logicOhmicNeg->SetVisAttributes(yellow);
  logicOhmicPos->SetVisAttributes(yellow);
  logicOptic->SetVisAttributes(grayc);


  if (mercuryGranularity)  logicGrain->SetVisAttributes(grayc);

  //always return the physical World
    
  PrintApparateParameters();

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorConstruction::ConstructSDandField()
{
   //                               
  // Sensitive Detectors 
  //
  if (HPGeSD.Get() == 0) 
    {    
      XrayFluoSD* SD = new XrayFluoSD ("HPGeSD",this);
      HPGeSD.Put( SD );
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(HPGeSD.Get());
  if (logicPixel)    
    SetSensitiveDetector(logicPixel,HPGeSD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The mercury is a sphere whose diamter is: "
	 << G4endl      
	 << mercuryDia/km
	 << " Km "
	 << G4endl
	 <<" Material: " << logicMercury->GetMaterial()->GetName() 
	 <<G4endl
	 <<"The Detector is a slice  " << DeviceThickness/(1.e-6*m) 
	 << " micron thick of " << pixelMaterial->GetName()<<G4endl
	 <<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorConstruction::UpdateGeometry()
{

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::Clean();
  G4LogicalVolumeStore::Clean();
  G4SolidStore::Clean();

  zRotPhiHPGe.rotateX(-1.*PhiHPGe);
  ComputeApparateParameters();  

  //Triggers a new call of Construct() and of all the geometry resets.
  G4RunManager::GetRunManager()->ReinitializeGeometry();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorConstruction::SetMercuryMaterial(G4String newMaterial)
{
  G4cout << "New Mercury Material: " << newMaterial << G4endl;
  logicMercury->SetMaterial(materials->GetMaterial(newMaterial));
  PrintApparateParameters();
   //GeometryHasBeenModified is called by the messenger
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











