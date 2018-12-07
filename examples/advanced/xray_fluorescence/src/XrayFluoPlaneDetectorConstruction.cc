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
// 29 aug 2003 Alfonso Mantero Created
// -------------------------------------------------------------------

#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorMessenger.hh"
#include "XrayFluoSD.hh"
#include "XrayFluoNistMaterials.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoPlaneDetectorConstruction::XrayFluoPlaneDetectorConstruction()
  : detectorType(0),planeGranularity(false), DeviceSizeX(0),
    DeviceSizeY(0),DeviceThickness(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidScreen(0),logicScreen(0),physiScreen(0),
    solidPlane (0),logicPlane(0),physiPlane (0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    screenMaterial(0),OhmicPosMaterial(0), OhmicNegMaterial(0),
    pixelMaterial(0),planeMaterial(0),
    defaultMaterial(0),HPGeSD(0)
  
{ 
  materials = XrayFluoNistMaterials::GetInstance();

  DefineDefaultMaterials();

  NbOfPixelRows     =  1; // should be 1
  NbOfPixelColumns  =  1; // should be 1
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
  PixelSizeXY       = 5 * cm; // should be 5
  PixelThickness = 3.5 * mm; //changed should be 3.5 mm

  G4cout << "PixelThickness(mm): "<< PixelThickness/mm << G4endl;
  G4cout << "PixelSizeXY(cm): "<< PixelSizeXY/cm << G4endl;

  ContactSizeXY  = 5 * cm; //should be the same as pixelSizeXY
  planeThickness = 5 * cm;
  planeSizeXY    = 5. * m;

  OhmicNegThickness = 0.005*mm;
  OhmicPosThickness = 0.005*mm;

  screenThickness = 5 * mm;

  ThetaHPGe = 0. * deg;
  PhiHPGe = 0. * deg;


  DistDe = 0.5 * m;

  distScreen = DistDe + (screenThickness+PixelThickness)/2+OhmicPosThickness ;

  grainDia = 1 * mm;


  PixelCopyNb=0;
  grainCopyNb=0;
  G4String defaultDetectorType = "sili";
  ComputeApparateParameters();
  SetDetectorType(defaultDetectorType);
  
  // create commands for interactive definition of the apparate
  
  detectorMessenger = new XrayFluoPlaneDetectorMessenger(this);
  G4cout << "XrayFluoPlaneDetectorConstruction created" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoPlaneDetectorConstruction* XrayFluoPlaneDetectorConstruction::instance = 0;

XrayFluoPlaneDetectorConstruction* XrayFluoPlaneDetectorConstruction::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoPlaneDetectorConstruction;
     
    }
  return instance;
}

void XrayFluoPlaneDetectorConstruction::SetDetectorType(G4String type)
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
      execp << type + "detector type unknown";
      G4Exception("XrayFluoPlaneDetectorConstruction::SetDetectorType()","example-xray_fluorescence03",
	  FatalException, execp);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoVDetectorType* XrayFluoPlaneDetectorConstruction::GetDetectorType() const 
{
  return detectorType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlaneDetectorConstruction::~XrayFluoPlaneDetectorConstruction()

{ 
  delete detectorMessenger;
  delete detectorType;
  G4cout << "XrayFluoPlaneDetectorConstruction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoPlaneDetectorConstruction::Construct()
{
  return ConstructApparate();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlaneDetectorConstruction::DefineDefaultMaterials()
{


  //define materials of the apparate

  planeMaterial = materials->GetMaterial("Anorthosite");
  screenMaterial = materials->GetMaterial("G4_Pb");
  pixelMaterial = materials->GetMaterial("G4_Si");
  OhmicPosMaterial = materials->GetMaterial("G4_Cu");
  OhmicNegMaterial = materials->GetMaterial("G4_Pb");
  defaultMaterial = materials->GetMaterial("G4_Galactic");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* XrayFluoPlaneDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition 
  
  //ComputeApparateParameters();
  
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

      z = -1. * DistDe; //* std::cos(ThetaHPGe);
      y = 0.*cm; //distScreen * std::sin(ThetaHPGe);
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
	z = DistDe * std::cos(ThetaHPGe);
	y =DistDe * std::sin(ThetaHPGe);
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
      z = -1 * distScreen; //* std::cos(ThetaHPGe);
      y = 0.*cm; //distScreen * std::sin(ThetaHPGe);
      x = 0.*cm;
      physiScreen = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)), 
				    "DetectorScreen",	//its name
				    logicScreen,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
    }

    //Plane

  if (planeGranularity) {

    solidPlane=0;  logicPlane=0;  physiPlane=0;
    if (planeThickness > 0.)  
      {
	solidPlane = new G4Box("Plane",		//its name
				planeSizeXY/2,planeSizeXY/2,planeThickness/2);//size
	
	logicPlane= new G4LogicalVolume(solidPlane,	//its solid
					 defaultMaterial,	//its material
					 "Plane");	//its name
	
	physiPlane = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Plane",	//its name
					logicPlane,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number
	
      }




    G4int nbOfGrainsX = ((G4int)(planeSizeXY/grainDia)) -1 ;
    
    // y dim of a max density plane is 2rn-(n-1)ar, wehere a = (1-(std::sqrt(3)/2)), n is 
    // number of rows and r the radius of the grain. so the Y-dim of the plane must 
    // be greater or equal to this. It results that nmust be <= (PlaneY-a)/(1-a).
    // Max Y shift of the planes superimposing along Z axis is minor (2/std::sqrt(3)r)

    G4double a = (1.-(std::sqrt(3.)/2.));
    G4int nbOfGrainsY =  (G4int) ( ((planeSizeXY/(grainDia/2.)) -a)/(2.-a) ) -1;

    // same for the z axis, but a = 2 * (std::sqrt(3) - std::sqrt(2))/std::sqrt(3)

    G4double b = 2. * (std::sqrt(3.) - std::sqrt(2.))/std::sqrt(3.);
    G4int nbOfGrainsZ =  (G4int) ( ((planeThickness/(grainDia/2.)) -b)/(2.-b) )-1;

    if (planeThickness > 0.){
      
      solidGrain=0; logicGrain=0; physiGrain=0;
      solidGrain = new G4Sphere("Grain",0.,			
				grainDia/2,0., twopi, 0., pi);
      
      logicGrain = new G4LogicalVolume(solidGrain,	
				       planeMaterial,	//its material
				       "Grain");	        //its name
      G4ThreeVector grainPosition; 
      G4double grainInitPositionX = 0;
      G4double grainInitPositionY = 0;
      G4double grainInitPositionZ = (-1.*planeThickness/2.+grainDia/2.);
      G4double grainStepX = grainDia = 0;
      G4double grainStepY = grainDia*(1.-(0.5-(std::sqrt(3.)/4.)));
      G4double grainStepZ = grainDia*std::sqrt(2./3.);
      
      for ( G4int k=0; k < nbOfGrainsZ ; k++ ) {
	for ( G4int j=0; j < nbOfGrainsY ; j++ ) {
	  for ( G4int i=0; i < nbOfGrainsX ; i++ ) {
	    
	    // Now we identify the layer and the row where the grain is , to place it in the right position
	    
	    
	    
	    if (k%3 == 0) { // first or (4-multiple)th layer: structure is ABCABC
	      grainInitPositionY = (-1.*planeSizeXY/2.+grainDia/2.);    
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*planeSizeXY/2.+grainDia/2.);
	      }
	      
	      else if ( ((j+1) % 2)  == 0 ) {
		grainInitPositionX = (-1.*planeSizeXY/2.+ grainDia);		
	      }
	      
	    }	      
	    else if ( ((k+2) % 3) == 0 ) { // B-layer
	      
	      grainInitPositionY = ( (-1.*planeSizeXY/2.) + (grainDia/2.)*(1. + (1./std::sqrt(3.)) ) );
	      
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*planeSizeXY/2.+grainDia);
	      }
	      
	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*planeSizeXY/2.+grainDia/2);		
	      }
	      
	    }
	    
	    else if ( (k+1)%3 == 0 ) { // B-layer
	      
	      grainInitPositionY = (-1.*planeSizeXY/2.+(grainDia/2.)*(1.+2./std::sqrt(3.)) );
	      
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*planeSizeXY/2.+grainDia/2.);
	      }
	      
	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*planeSizeXY/2.+grainDia);		
	      }
	      
	    }
	    
	    physiGrain = new G4PVPlacement(0,	       
					   G4ThreeVector( grainInitPositionX + i*grainStepX, 
							  grainInitPositionY + j*grainStepY,
							  grainInitPositionZ + k*grainStepZ),
					   "Grain",  
					   logicGrain,	 //its logical volume
					   physiPlane, //its mother  volume
					   false,	 //no boolean operation
					   grainCopyNb);//copy number    
	    
	    grainCopyNb = grainCopyNb +1; 
	  }
	}
      }
    }    
  }
  else {     
      
    solidPlane=0;  logicPlane=0;  physiPlane=0;
    if (planeThickness > 0.)  
      {
	solidPlane = new G4Box("Plane",		//its name
				planeSizeXY/2,planeSizeXY/2,planeThickness/2);//size
	  
	logicPlane= new G4LogicalVolume(solidPlane,	//its solid
					 planeMaterial,	//its material
					 "Plane");	//its name
	  
	physiPlane = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Plane",	//its name
					logicPlane,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number
	  
      }  
  }

  // Visualization attributes
  
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());
   G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255. ,255/255. ,51/255. ));
   G4VisAttributes * red= new G4VisAttributes( G4Colour(255/255. , 0/255. , 0/255. ));
   G4VisAttributes * blue= new G4VisAttributes( G4Colour(0/255. , 0/255. ,  255/255. ));
   G4VisAttributes * grayc= new G4VisAttributes( G4Colour(128/255. , 128/255. ,  128/255. ));
   G4VisAttributes * lightGray= new G4VisAttributes( G4Colour(178/255. , 178/255. ,  178/255. ));
  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  red->SetVisibility(true);
  red->SetForceSolid(true);
  blue->SetVisibility(true);
  grayc->SetVisibility(true);
  grayc->SetForceSolid(true);
  lightGray->SetVisibility(true);
  lightGray->SetForceSolid(true);
  simpleBoxVisAtt->SetVisibility(true);
 
  logicPixel->SetVisAttributes(red); //modified!!!
  logicHPGe->SetVisAttributes(blue);

  logicPlane->SetVisAttributes(lightGray);
  

  logicScreen->SetVisAttributes(grayc);
  logicOhmicNeg->SetVisAttributes(yellow);
  logicOhmicPos->SetVisAttributes(yellow);



  if (planeGranularity)  logicGrain->SetVisAttributes(grayc);

  //always return the physical World
    
  PrintApparateParameters();

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...

void XrayFluoPlaneDetectorConstruction::ConstructSDandField()
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

void XrayFluoPlaneDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The plane is a box whose size is: "
	 << G4endl      
	 << planeThickness/cm
	 << " cm * "
	 << planeSizeXY/cm
	 << " cm * "
	 << planeSizeXY/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << logicPlane->GetMaterial()->GetName() 
	 <<G4endl
	  <<"The Detector is a slice  " << DeviceThickness/(1.e-6*m) <<  " micron thick of " << pixelMaterial->GetName()
	 <<G4endl
	 

<<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlaneDetectorConstruction::UpdateGeometry()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::Clean();
  G4LogicalVolumeStore::Clean();
  G4SolidStore::Clean();
 
  zRotPhiHPGe.rotateX(-1.*PhiHPGe);
  //Triggers a new call of Construct() and of all the geometry resets.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlaneDetectorConstruction::DeleteGrainObjects()
{
  if (planeGranularity) { 
    delete solidGrain; 
    delete logicGrain;
    delete physiGrain;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlaneDetectorConstruction::SetPlaneMaterial(G4String newMaterial)
{
  //G4cout << "Material!!!!" << newMaterial << G4endl;
  logicPlane->SetMaterial(materials->GetMaterial(newMaterial));
  PrintApparateParameters();
  //GeometryHasBeenModified is called by the messenger
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....












