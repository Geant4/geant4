//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "myDetectorConstruction.hh"
#include "myDetectorMessenger.hh"
#include "mySiSD.hh"
#include "myHPGeSD.hh"
#include "mySampleSD.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

myDetectorConstruction::myDetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidSi(NULL),logicSi(NULL),physiSi(NULL),
 solidHPGe(NULL),logicHPGe(NULL),physiHPGe(NULL),
 solidSample (NULL),logicSample(NULL),physiSample (NULL),
 solidDia1(NULL),logicDia1(NULL),physiDia1(NULL),
 solidDia2(NULL),logicDia2(NULL),physiDia2(NULL),
 solidDia3(NULL),logicDia3(NULL),physiDia3(NULL),
 SiMaterial(NULL),HPGeMaterial(NULL),sampleMaterial(NULL),
 Dia1Material(NULL),Dia2Material(NULL),Dia3Material(NULL),
 defaultMaterial(NULL),
 sampleSD(NULL),SiSD(NULL),HPGeSD(NULL)
{
  
  SiSizeYZ = 1. * cm;
  SiThickness = 1.2 * cm;
  HPGeSizeYZ = 1. * cm;
  HPGeThickness = 1.2 * cm;
  SampleThickness = 1. * mm;
  SampleSizeYZ = 3. * cm;
  Dia1Thickness = 1. *mm;
  Dia2Thickness = 1. *mm;
  Dia3Thickness = 1. *mm;
  Dia1SizeYZ = 1. *cm;
  Dia2SizeYZ = 1. *cm;
  Dia3SizeYZ = 1. *cm;
  DiaInnerSize = 0.5 * cm;
 
  ThetaHPGe = 150. * deg;
  ThetaSi = 210. * deg;
  DistDe = 77.0 * mm;
  PhiHPGe = 150. * deg;
  PhiSi = 210. * deg;
  ThetaDia1 = 150. * deg;
  ThetaDia2 = 210. * deg;
  ThetaDia3 = 180. * deg;
  PhiDia3 = 90. * deg;
  DistDia = 67.0 * mm;
  PhiDia1 = 90. * deg;
  PhiDia2 = 90. * deg;
  AlphaDia1 = 150. * deg;
  AlphaDia2 = 210. * deg;
  AlphaDia3 = 180. * deg;
  ComputeApparateParameters();
 
  // create commands for interactive definition of the apparate

  detectorMessenger = new myDetectorMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

myDetectorConstruction::~myDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* myDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructApparate();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::DefineMaterials()
{
  //define elements
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
  
 
 a = 74.9216 * g/mole;
  G4Element * As = new G4Element( name="arsenic",symbol="As",z= 33.,a);
 
  a = 69.72* g/mole;
  G4Element * Ga = new G4Element(name="gallium",symbol="Ga",z= 31.,a);

  a = 55.85*g/mole;
  G4Element* Fe = new G4Element(name="Iron"  ,symbol="Fe", z=26., a);

  density = 7.86 * g/cm3;
  
  G4int  natoms,ncomponents;
  G4double temperature, pressure;
  
  G4Material * FeMaterial = new G4Material(name="Iron",density,ncomponents=1);
  FeMaterial->AddElement(Fe,natoms=1);
 
  //define gallium arsenide

 density = 5.32 * g/cm3;
 G4Material * GaAs = new G4Material(name ="gallium arsenide",density,ncomponents=2);
 GaAs->AddElement(Ga,natoms=1);
 GaAs->AddElement(As,natoms=1);
 
 //define silicon

  density = 2.333*g/cm3;
  a = 28.09*g/mole;
  G4Material* Si = new G4Material(name="Silicon",z=14., a,density);
  
  //definr carbon
  
  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  //define hydrogen

 a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
   
//define scintillator

 density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);

  //define aluminium

density = 2.700*g/cm3;
a = 26.98*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

//define vacuum

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material * Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the apparate

 SiMaterial = Si;
 sampleMaterial = FeMaterial;
 defaultMaterial = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* myDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition 
  
ComputeApparateParameters();

//world

solidWorld = new G4Box("World",	      		        //its name
                 WorldSizeX/2,WorldSizeYZ/2,WorldSizeYZ/2);	//its size

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

//SiDetector

 solidSi = NULL;  physiSi = NULL;  logicSi=NULL;

if (SiThickness > 0.)  
{
   solidSi = new G4Box("SiDetector",		//its name
  		       SiThickness/2,SiSizeYZ/2,SiSizeYZ/2);//size

   
   logicSi = new G4LogicalVolume(solidSi,	//its solid
      				       SiMaterial,	//its material
      				       "SiDetector");	//its name
  		      
   zRotPhiSi.rotateZ(PhiSi);
   G4double x,y,z;
   x = DistDe * cos(ThetaSi);
   y =DistDe * sin(ThetaSi);
   z = 0.*cm;
   physiSi = new G4PVPlacement(G4Transform3D(zRotPhiSi,G4ThreeVector(x,y,z)),                                           "SiDetector",	//its name
                                   logicSi,	//its logical volume
                                   physiWorld,	//its mother  volume
                                   false,		//no boolean operation
                                   0);		//copy number
}
//HPGeDetector

 solidHPGe = NULL;  physiHPGe = NULL;  logicHPGe=NULL;

if (HPGeThickness > 0.)  
{
   solidHPGe = new G4Box("HPGeDetector",		//its name
  		       HPGeThickness/2,HPGeSizeYZ/2,HPGeSizeYZ/2);//size

   
   logicHPGe = new G4LogicalVolume(solidHPGe,	//its solid
      				       HPGeMaterial,	//its material
      				       "HPGeDetector");	//its name
  		      
   zRotPhiHPGe.rotateZ(PhiHPGe);
   G4double x,y,z;
   x = DistDe * cos(ThetaHPGe);
   y =DistDe * sin(ThetaHPGe);
   z = 0.*cm;
   physiHPGe = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)),                                           "HPGeDetector",	//its name
                                   logicHPGe,	//its logical volume
                                   physiWorld,	//its mother  volume
                                   false,		//no boolean operation
                                   0);		//copy number
}     
 
//Sample

  solidSample=NULL;  logicSample=NULL;  physiSample=NULL;

if (SampleThickness > 0.)  
{
   solidSample = new G4Box("Sample",		//its name
  		       SampleThickness/2,SampleSizeYZ/2,SampleSizeYZ/2);//size
   
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
 //Diaphragm1

 solidDia1 = NULL;  physiDia1 = NULL;  logicDia1=NULL;

if (Dia1Thickness > 0.)  
{
   solidDia1 = new G4Tubs("Diaphragm1",		//its name
  		       DiaInnerSize/2,
			  Dia1SizeYZ/2,
			  Dia1Thickness/2,
			  0,
			  360);//size

   
   logicDia1 = new G4LogicalVolume(solidDia1,	//its solid
      				       Dia1Material,	//its material
      				       "Diaphragm1");	//its name
  		      
   zRotPhiDia1.rotateY(PhiDia1);
   zRotPhiDia1.rotateZ(AlphaDia1);
   G4double x,y,z;
   x = DistDia * cos(ThetaDia1);
   y =DistDia * sin(ThetaDia1);
   z = 0.*cm;
   physiDia1 = new G4PVPlacement(G4Transform3D(zRotPhiDia1,G4ThreeVector(x,y,z)),                                           "Diaphragm1",	//its name
                                   logicDia1,	//its logical volume
                                   physiWorld,	//its mother  volume
                                   false,		//no boolean operation
                                   0);		//copy number
}     
 //Diaphragm2

 solidDia2 = NULL;  physiDia2 = NULL;  logicDia2=NULL;

if (Dia2Thickness > 0.)  
{
   solidDia2 = new G4Tubs("Diaphragm2",
			  DiaInnerSize/2,
			  Dia2SizeYZ/2,
			  Dia2Thickness/2,
			  0,
			  360);

   
   logicDia2 = new G4LogicalVolume(solidDia2,	//its solid
      				       Dia2Material,	//its material
      				       "Diaphragm2");	//its name
  		      
   zRotPhiDia2.rotateY(PhiDia2);
   zRotPhiDia2.rotateZ(AlphaDia2);
   G4double x,y,z;
   x = DistDia * cos(ThetaDia2);
   y =DistDia * sin(ThetaDia2);
   z = 0.*cm;
   physiDia2 = new G4PVPlacement(G4Transform3D(zRotPhiDia2,G4ThreeVector(x,y,z)),                                           "Diaphragm1",	//its name
                                   logicDia2,	//its logical volume
                                   physiWorld,	//its mother  volume
                                   false,		//no boolean operation
                                   0);		//copy number
} 
//Diaphragm3

 solidDia3 = NULL;  physiDia3 = NULL;  logicDia3 =NULL;

if (Dia3Thickness > 0.)  
{
   solidDia3 = new G4Tubs("Diaphragm3",
			 DiaInnerSize/2,
			  Dia3SizeYZ/2,
			  Dia3Thickness/2,
			  0,
			  360);

   
   logicDia3 = new G4LogicalVolume(solidDia3,	//its solid
      				       Dia3Material,	//its material
      				       "Diaphragm2");	//its name
  		      
   zRotPhiDia3.rotateY(PhiDia3);
 zRotPhiDia3.rotateZ(AlphaDia3);
   G4double x,y,z;
   x = DistDia * cos(ThetaDia3);
   y =DistDia * sin(ThetaDia3);
   z = 0.*cm;
   physiDia3 = new G4PVPlacement(G4Transform3D(zRotPhiDia3,G4ThreeVector(x,y,z)),                                           "Diaphragm1",	//its name
                                   logicDia3,	//its logical volume
                                   physiWorld,	//its mother  volume
                                   false,		//no boolean operation
                                   0);		//copy number
}    
    
G4SDManager* SDman = G4SDManager::GetSDMpointer();

if(!SiSD)
{
  SiSD = new mySiSD ("SiSD",this);
  SDman->AddNewDetector(SiSD);
}

if(!sampleSD)
{
sampleSD = new mySampleSD ("SampleSD",this);
  SDman->AddNewDetector(sampleSD);
}

if(!HPGeSD)
{
  HPGeSD = new myHPGeSD ("HPGeSD",this);
  SDman->AddNewDetector(HPGeSD);
}
if (logicSi)
{
logicSi->SetSensitiveDetector(SiSD);
}

if (logicHPGe)
{
logicHPGe->SetSensitiveDetector(HPGeSD);
}

if (logicSample)
{
    logicSample->SetSensitiveDetector(sampleSD); 
}
// Visualization attributes

  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicSi->SetVisAttributes(simpleBoxVisAtt);
  logicHPGe->SetVisAttributes(simpleBoxVisAtt);
  logicSample->SetVisAttributes(simpleBoxVisAtt);
  logicDia1->SetVisAttributes(simpleBoxVisAtt);
  logicDia2->SetVisAttributes(simpleBoxVisAtt);
  logicDia3->SetVisAttributes(simpleBoxVisAtt);

  //always return the physical World

PrintApparateParameters();
return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The sample is a box whose size is: "
	 << G4endl      
	 << SampleThickness/cm
	 << " cm * "
	 << SampleSizeYZ/cm
	 << " cm * "
	 << SampleSizeYZ/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << sampleMaterial->GetName() 
	 <<G4endl
	 <<"The SiDetector is a slice  " << SiThickness/(1.e-6*m) <<  " micron thick"
	 <<G4endl
	 <<"Size: "<< SiSizeYZ/cm<< " cm * "<<SiSizeYZ/cm 
	 <<G4endl
	  <<"The HPGeDetector is a slice  " << HPGeThickness/(1.e-6*m) <<  " micron thick"
	 <<G4endl
	 <<"Size: "<< HPGeSizeYZ/cm<< " cm * "<<HPGeSizeYZ/cm 
	 <<G4endl

<<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSampleMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {sampleMaterial = pttoMaterial;
      logicSample->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSiMaterial(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {SiMaterial = pttoMaterial;
      logicSi->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetHPGeMaterial(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {HPGeMaterial = pttoMaterial;
      logicHPGe->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia1Material(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {Dia1Material = pttoMaterial;
      logicDia1->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia2Material(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {Dia2Material = pttoMaterial;
      logicDia2->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia3Material(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {Dia3Material = pttoMaterial;
      logicDia3->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSampleThickness(G4double val)
{
  // change Sample thickness 
   SampleThickness = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSiThickness(G4double val)
{
  // change Sensor thickness 
  SiThickness = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetHPGeThickness(G4double val)
{
  // change Sensor thickness 
  HPGeThickness = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia1Thickness(G4double val)
{
  // change Sensor thickness 
  Dia1Thickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia2Thickness(G4double val)
{
  // change Sensor thickness 
  Dia2Thickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia3Thickness(G4double val)
{
  // change Sensor thickness 
  Dia3Thickness = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSampleSizeYZ(G4double val)
{
 // change the transverse size of the sample and recompute the world size
  SampleSizeYZ = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSiSizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  SiSizeYZ = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetHPGeSizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  HPGeSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia1SizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  Dia1SizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia2SizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  Dia2SizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia3SizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  Dia3SizeYZ = val;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSensorDistance(G4double val)
{
  // change the distance between the sensor and the sample
  DistDe = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDiaphragmDistance(G4double val)
{
  // change the distance between the sensor and the sample
  DistDia = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSiAzimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  ThetaSi = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetHPGeAzimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  ThetaHPGe = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia1Azimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  ThetaDia1 = val;
} 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia2Azimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  ThetaDia2 = val;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia3Azimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  ThetaDia3 = val;
} 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetSiRotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  PhiSi = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetHPGeRotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  PhiHPGe = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia1Rotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  PhiDia1 = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void myDetectorConstruction::SetDia2Rotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  PhiDia2 = val;
}

void myDetectorConstruction::SetDia3Rotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  PhiDia3 = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void myDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructApparate());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
