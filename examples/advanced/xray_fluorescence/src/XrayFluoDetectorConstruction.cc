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
// $Id: XrayFluoDetectorConstruction.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//    Nov 2002 Alfonso Mantero materials added, 
//             Material selection implementation
// 16 Jul 2003 Alfonso Mantero Detector type selection added + minor fixes
// -------------------------------------------------------------------

#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoSD.hh"
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
#include "XrayFluoNistMaterials.hh"
#include "G4SDManager.hh"

// #include "G4Region.hh"
// #include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoDetectorConstruction::XrayFluoDetectorConstruction()
  : aNavigator(0), detectorType(0),sampleGranularity(false), phaseSpaceFlag(false),
    DeviceSizeX(0), DeviceSizeY(0),DeviceThickness(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidSample (0),logicSample(0),physiSample (0),
    solidDia1(0),logicDia1(0),physiDia1(0),
    solidDia3(0),logicDia3(0),physiDia3(0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidWindow(0), logicWindow(0), physiWindow(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    OhmicPosMaterial(0), OhmicNegMaterial(0),
    pixelMaterial(0),sampleMaterial(0),
    Dia1Material(0),Dia3Material(0),
    defaultMaterial(0), windowMaterial (0)  
{ 
  materials = XrayFluoNistMaterials::GetInstance();

  HPGeSD.Put(0);

  aNavigator = new G4Navigator();
 
  DefineDefaultMaterials();

  NbOfPixelRows     =  1; // should be 1
  NbOfPixelColumns  =  1; // should be 1
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
  PixelSizeXY       =  std::sqrt(40.) * mm;
  PixelThickness = 2.7 * mm; //should be 3.5 mm

  G4cout << "PixelThickness(mm): "<< PixelThickness/mm << G4endl;
  G4cout << "PixelSizeXY(cm): "<< PixelSizeXY/cm << G4endl;

  ContactSizeXY     = PixelSizeXY; //std::sqrt(40) * mm; //should be the same as PixelSizeXY
  SampleThickness = 4 * mm;
  SampleSizeXY = 3. * cm;
  Dia1Thickness = 1. *mm;
  Dia3Thickness = 1. *mm;
  Dia1SizeXY = 3. *cm;
  Dia3SizeXY = 3. *cm;


  DiaInnerSize = 2.9 * cm; //(Hole in the detector's diaphragm) it was 1 mm


  OhmicNegThickness = 1e-6*cm;// 0.005
  OhmicPosThickness = 1e-6*cm;// 0.005
  windowThickness = 0.008 * cm; //value for aif detector
  ThetaHPGe = 135. * deg;
  PhiHPGe = 225. * deg;

  ThetaDia1 = 135. * deg;
  PhiDia1 = 90. * deg;
  AlphaDia1 = 225. * deg;

  AlphaDia3 = 180. * deg;
  Dia3Dist =  66.5 * mm;
  Dia3InnerSize = 1. * mm;
  ThetaDia3 = 180. * deg;
  PhiDia3 = 90. * deg;

  DistDia = 66.5 * mm;
  DistDe =DistDia+ (Dia1Thickness
		    +PixelThickness)/2+OhmicPosThickness+windowThickness ;

  grainDia = 1 * mm;
  PixelCopyNb=0;
  grainCopyNb=0;
  G4String defaultDetectorType = "sili";
  ComputeApparateParameters();

//   G4String regName = "SampleRegion";
//   sampleRegion = new G4Region(regName);  

  if (!phaseSpaceFlag) SetDetectorType(defaultDetectorType);
  
  // create commands for interactive definition of the apparate
  
  detectorMessenger = new XrayFluoDetectorMessenger(this);

  G4cout << "XrayFluoDetectorConstruction created" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


XrayFluoDetectorConstruction* XrayFluoDetectorConstruction::instance = 0;

XrayFluoDetectorConstruction* XrayFluoDetectorConstruction::GetInstance()
{
  if (instance == 0)
    {
      instance = new XrayFluoDetectorConstruction;
     
    }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetDetectorType(G4String type) 
{
  if (type=="sili")
    {
      detectorType = XrayFluoSiLiDetectorType::GetInstance();
    }
   else if (type=="hpge")
     {
       detectorType = XrayFluoHPGeDetectorType::GetInstance();
    }/*
   else if (type=="aifira")
     {
       detectorType = XrayFluoAifSiLi::GetInstance();
       }*/
  else 
    {
      G4ExceptionDescription execp;
      execp << type + "detector type unknown";
      G4Exception("XrayFluoDataSet::LoadData()","example-xray_fluorescence06",
	  FatalException, execp);
    }
  //GeometryHasBeenModified invoked by the messenger
  
}

XrayFluoVDetectorType* XrayFluoDetectorConstruction::GetDetectorType() const
{
  return detectorType;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::~XrayFluoDetectorConstruction()

{ 
  delete detectorMessenger;
  delete detectorType;
  G4cout << "XrayFluoDetectorConstruction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoDetectorConstruction::Construct()
{
  return ConstructApparate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DefineDefaultMaterials()
{


  //define materials of the apparate

  sampleMaterial = materials->GetMaterial("Dolorite");
  Dia1Material = materials->GetMaterial("G4_Pb");
  Dia3Material = materials->GetMaterial("G4_Galactic");
  pixelMaterial = materials->GetMaterial("SiLi");
  //OhmicPosMaterial = materials->GetMaterial("G4_Cu");
  OhmicPosMaterial = materials->GetMaterial("G4_Ni");
  OhmicNegMaterial = materials->GetMaterial("G4_Pb");
  defaultMaterial = materials->GetMaterial("G4_Galactic");
  windowMaterial = materials->GetMaterial("G4_Be");
}

void XrayFluoDetectorConstruction::SetOhmicPosThickness(G4double val)
{
  
  if (!phaseSpaceFlag) {    
    
    
    if (val == 0.0) {
      OhmicPosMaterial = materials->GetMaterial("G4_Galactic");
    }
    else {
      OhmicPosThickness = val;
      //OhmicPosMaterial = materials->GetMaterial("G4_Cu");
      OhmicPosMaterial = materials->GetMaterial("G4_Ni");
    }
    
  }
  else{
    G4cout << "Not available in this configuration" << G4endl;
  }
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* XrayFluoDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition 
  
  //ComputeApparateParameters();

  //world and associated navigator
  
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

  aNavigator->SetWorldVolume(physiWorld);

 
  //HPGeDetector

  if (!phaseSpaceFlag) {
    
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
	z = DistDe * std::cos(ThetaHPGe);
	y =DistDe * std::sin(ThetaHPGe);
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

	  /////////// widow place here! ////////////////
	  // OhmicPos
	  solidWindow=0; logicWindow=0; physiWindow=0;  
	  
	  if (windowThickness > 0.) 
	    { solidWindow = new G4Box("Window",		//its name
					PixelSizeXY/2,PixelSizeXY/2,windowThickness/2); 
	    
	    logicWindow = new G4LogicalVolume(solidWindow,    //its solid
						windowMaterial, //its material
						"Window");      //its name
	    
	    physiWindow = new G4PVPlacement(0,	
					      G4ThreeVector(0.,
							    0.,
							    ((-PixelThickness-windowThickness)/2)
							    -OhmicPosThickness),  
					      "OhmicWindow",  
					      logicWindow,
					      physiHPGe,  
					      false,     
					      PixelCopyNb); 
	    
	    }


	
	  PixelCopyNb += PixelCopyNb; 
	  G4cout << "PixelCopyNb: " << PixelCopyNb << G4endl;
	}
      
      }
    
  }
  
  //Sample
  
  if (sampleGranularity) {

    solidSample=0;  logicSample=0;  physiSample=0;
    if (SampleThickness > 0.)  
      {
	solidSample = new G4Box("Sample",		//its name
				SampleSizeXY/2,SampleSizeXY/2,SampleThickness/2);//size
	
	logicSample= new G4LogicalVolume(solidSample,	//its solid
					 defaultMaterial,	//its material
					 "Sample");	//its name
	
	physiSample = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Sample",	//its name
					logicSample,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number
	
      }




    G4int nbOfGrainsX = ((G4int)(SampleSizeXY/grainDia)) -1 ;
    
    // y dim of a max density plane is 2rn-(n-1)ar, wehere a = (1-(std::sqrt(3)/2)), n is 
    // number of rows and r the radius of the grain. so the Y-dim of the sample must 
    // be greater or equal to this. It results that nmust be <= (SampleY-a)/(1-a).
    // Max Y shift of the planes superimposing along Z axis is minor (2/std::sqrt(3)r)

    G4double a = (1.-(std::sqrt(3.)/2.));
    G4int nbOfGrainsY =  (G4int) ( ((SampleSizeXY/(grainDia/2.)) -a)/(2.-a) ) -1;

    // same for the z axis, but a = 2 * (std::sqrt(3) - std::sqrt(2))/std::sqrt(3)

    G4double b = 2. * (std::sqrt(3.) - std::sqrt(2.))/std::sqrt(3.);
    G4int nbOfGrainsZ =  (G4int) ( ((SampleThickness/(grainDia/2.)) -b)/(2.-b) )-1;

    if (SampleThickness > 0.){
      
      solidGrain=0; logicGrain=0; physiGrain=0;
      solidGrain = new G4Sphere("Grain",0.,			
				grainDia/2,0., twopi, 0., pi);
      
      logicGrain = new G4LogicalVolume(solidGrain,	
				       sampleMaterial,	//its material
				       "Grain");	        //its name
      G4ThreeVector grainPosition; 
      G4double grainInitPositionX = 0;
      G4double grainInitPositionY = 0;
      G4double grainInitPositionZ = (-1.*SampleThickness/2.+grainDia/2.);
      G4double grainStepX = grainDia = 0;
      G4double grainStepY = grainDia*(1.-(0.5-(std::sqrt(3.)/4.)));
      G4double grainStepZ = grainDia*std::sqrt(2./3.);
      
      for ( G4int k=0; k < nbOfGrainsZ ; k++ ) {
	for ( G4int j=0; j < nbOfGrainsY ; j++ ) {
	  for ( G4int i=0; i < nbOfGrainsX ; i++ ) {
	    
	    // Now we identify the layer and the row where the grain is , to place it in the right position
	    
	    
	    
	    if (k%3 == 0) { // first or (4-multiple)th layer: structure is ABCABC
	      grainInitPositionY = (-1.*SampleSizeXY/2.+grainDia/2.);    
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2.);
	      }
	      
	      else if ( ((j+1) % 2)  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+ grainDia);		
	      }
	      
	    }	      
	    else if ( ((k+2) % 3) == 0 ) { // B-layer
	      
	      grainInitPositionY = ( (-1.*SampleSizeXY/2.) + (grainDia/2.)*(1. + (1./std::sqrt(3.)) ) );
	      
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia);
	      }
	      
	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2);		
	      }
	      
	    }
	    
	    else if ( (k+1)%3 == 0 ) { // B-layer
	      
	      grainInitPositionY = (-1.*SampleSizeXY/2.+(grainDia/2.)*(1.+2./std::sqrt(3.)) );
	      
	      if (j%2 ==0) { //first or (3-multiple)th row
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia/2.);
	      }
	      
	      else if ( (j+1)%2  == 0 ) {
		grainInitPositionX = (-1.*SampleSizeXY/2.+grainDia);		
	      }
	      
	    }
	    
	    physiGrain = new G4PVPlacement(0,	       
					   G4ThreeVector( grainInitPositionX + i*grainStepX, 
							  grainInitPositionY + j*grainStepY,
							  grainInitPositionZ + k*grainStepZ),
					   "Grain",  
					   logicGrain,	 //its logical volume
					   physiSample, //its mother  volume
					   false,	 //no boolean operation
					   grainCopyNb);//copy number    
	    
	    grainCopyNb = grainCopyNb +1; 
	  }
	}
      }
    }    
  }
  else {     
      
    solidSample=0;  logicSample=0;  physiSample=0;
    if (SampleThickness > 0.)  
      {
	solidSample = new G4Box("Sample",		//its name
				SampleSizeXY/2,SampleSizeXY/2,SampleThickness/2);//size
	  
	logicSample= new G4LogicalVolume(solidSample,	//its solid
					 sampleMaterial,	//its material
					 "Sample");	//its name
	  
	physiSample = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Sample",	//its name
					logicSample,	//its logical volume
					physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number
	  
      }  
  }

  if (!phaseSpaceFlag) {    
    //Diaphragm1
    
    solidDia1 = 0;  physiDia1 = 0;  logicDia1=0;
    
    if (Dia1Thickness > 0.)  
      {
	solidDia1 = new G4Tubs("Diaphragm1",		//its name
			       DiaInnerSize/2,
			       Dia1SizeXY/2,
			       Dia1Thickness/2,
			       0,
			       360*deg);//size
	
	
	logicDia1 = new G4LogicalVolume(solidDia1,	//its solid
					Dia1Material,	//its material
					"Diaphragm1");	//its name
	
	zRotPhiDia1.rotateX(AlphaDia1);
	G4double x,y,z;
	z = DistDia * std::cos(ThetaDia1);
	y =DistDia * std::sin(ThetaDia1);
	x = 0.*cm;
	physiDia1 = new G4PVPlacement(G4Transform3D(zRotPhiDia1,G4ThreeVector(x,y,z)),
				      "Diaphragm1",	//its name
				      logicDia1,	//its logical volume
				      physiWorld,	//its mother  volume
				      false,		//no boolean operation
				      0);		//copy number
      }  
    
    //Diaphragm3
    
    solidDia3 = 0;  physiDia3 = 0;  logicDia3 =0;
    
    if (Dia3Thickness > 0.)  
      {
	solidDia3 = new G4Tubs("Diaphragm3",
			       Dia3InnerSize/2,
			       Dia3SizeXY/2,
			       Dia3Thickness/2,
			       0,
			       360*deg);
	
	
      logicDia3 = new G4LogicalVolume(solidDia3,	//its solid
				      Dia3Material,	//its material
				      "Diaphragm3");	//its name
      
      zRotPhiDia3.rotateX(AlphaDia3);
      G4double x,y,z;
      z = Dia3Dist * std::cos(ThetaDia3);
      y =Dia3Dist * std::sin(ThetaDia3);
      x = 0.*cm;
      physiDia3 = new G4PVPlacement(G4Transform3D(zRotPhiDia3,G4ThreeVector(x,y,z)),                                           "Diaphragm3",	//its name
				    logicDia3,	//its logical volume
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
  G4VisAttributes * green= new G4VisAttributes( G4Colour(0/255. , 255/255. ,  0/255. ));

  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  red->SetVisibility(true);
  red->SetForceSolid(true);
  blue->SetVisibility(true);
  green->SetVisibility(true);
  green->SetForceSolid(true);
  grayc->SetVisibility(true);
  grayc->SetForceSolid(true);
  lightGray->SetVisibility(true);
  lightGray->SetForceSolid(true);
  simpleBoxVisAtt->SetVisibility(true);
  if (!phaseSpaceFlag) {  
    logicPixel->SetVisAttributes(red); //modified!!!
    logicHPGe->SetVisAttributes(blue);
    
    logicDia1->SetVisAttributes(lightGray);
    logicDia3->SetVisAttributes(lightGray);
    
    logicOhmicNeg->SetVisAttributes(yellow);
    logicOhmicPos->SetVisAttributes(yellow);
    
    logicWindow->SetVisAttributes(green);

  }
  logicSample->SetVisAttributes(simpleBoxVisAtt);

  if (sampleGranularity) logicSample->SetVisAttributes(simpleBoxVisAtt); // mandatory



  if (sampleGranularity)  logicGrain->SetVisAttributes(grayc);

  //always return the physical World
    
  PrintApparateParameters();

  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::ConstructSDandField()
{
  if (!phaseSpaceFlag) 
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
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The sample is a box whose size is: "
	 << G4endl      
	 << SampleThickness/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << logicSample->GetMaterial()->GetName() 
	 <<G4endl;
  if (!phaseSpaceFlag) {  
    G4cout <<"The Detector is a slice  " << DeviceThickness/(1.e-6*m) <<  " micron thick of " << pixelMaterial->GetName()
	   <<G4endl
	   << "The Anode is a slice " << OhmicPosThickness/mm << "mm thick of "<< OhmicPosMaterial->GetName()
	   <<G4endl;
      }
  G4cout <<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::UpdateGeometry()
{
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::Clean();
  G4LogicalVolumeStore::Clean();
  G4SolidStore::Clean();
 
  if (sampleRegion) 
    sampleRegion->RemoveRootLogicalVolume(logicSample);
   
  zRotPhiHPGe.rotateX(-1.*PhiHPGe);
  zRotPhiDia1.rotateX(-1.*AlphaDia1);
  zRotPhiDia3.rotateX(-1.*AlphaDia3);

  //Triggers a new call of Construct() and of all the geometry resets.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DeleteGrainObjects()
{
  if (sampleGranularity) { 
    delete solidGrain; 
    delete logicGrain;
    delete physiGrain;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector XrayFluoDetectorConstruction::GetDetectorPosition() const
{
  G4double 	z = DistDe * std::cos(ThetaHPGe);
  G4double	y = DistDe * std::sin(ThetaHPGe);
  G4double	x = 0.*cm;

  G4ThreeVector position(x,y,z);

  return position;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSampleMaterial(G4String newMaterial)
{
    G4cout << "Material Change in Progress " << newMaterial << G4endl;
    sampleMaterial = materials->GetMaterial(newMaterial);
    logicSample->SetMaterial(sampleMaterial);
    PrintApparateParameters();
    //GeometryHasBeenModified is called by the messenger
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











