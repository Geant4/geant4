//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoSensorSD.hh"
#include "XrayFluoSampleSD.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::XrayFluoDetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidSensor(NULL),logicSensor(NULL),physiSensor(NULL),
 solidSample (NULL),logicSample(NULL),physiSample (NULL),
 sensorMaterial(NULL),sampleMaterial(NULL),defaultMaterial(NULL),
 sensorSD(NULL),sampleSD(NULL)
{
  //default parameter values (sizes)
  //SensorThickness = 0.25 * mm;
  SensorThickness = 10. * cm;
  //SensorSizeYZ = 1.5 * mm;
  SensorSizeYZ = 10. * cm;
  SampleThickness = 5. * cm;
  SampleSizeYZ = 2. * cm;
 
//default parameter values (position of the sensor)
  Theta = 30. * deg;
  Dist = 0.5 * m;
  Phi = 30. * deg;
  ComputeApparateParameters();
 
  // create commands for interactive definition of the apparate

  detectorMessenger = new XrayFluoDetectorMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::~XrayFluoDetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructApparate();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DefineMaterials()
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
//define vacuum

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material * Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

G4cout << *(G4Material::GetMaterialTable()) << G4endl;

//default materials of the apparate

 sensorMaterial = GaAs;
 sampleMaterial = FeMaterial;
 defaultMaterial = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* XrayFluoDetectorConstruction::ConstructApparate()
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

//sensor

 solidSensor = NULL;  physiSensor = NULL;  logicSensor=NULL;

if (SensorThickness > 0.)  
{
   solidSensor = new G4Box("Sensor",		//its name
  		       SensorThickness/2,SensorSizeYZ/2,SensorSizeYZ/2);//size

   
   logicSensor = new G4LogicalVolume(solidSensor,	//its solid
      				       sensorMaterial,	//its material
      				       "Sensor");	//its name
  		      
   zRotPhi.rotateZ(Phi);
   G4double x,y,z;
   x = Dist * cos(Theta);
   y =Dist * sin(Theta);
   z = 0.*cm;
   physiSensor = new G4PVPlacement(G4Transform3D(zRotPhi,G4ThreeVector(x,y,z)),                                           "Sensor",	//its name
                                   logicSensor,	//its logical volume
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
   
   logicSample = new G4LogicalVolume(solidSample,	//its solid
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
G4SDManager* SDman = G4SDManager::GetSDMpointer();

if(!sensorSD)
{
  sensorSD = new XrayFluoSensorSD ("SensorSD",this);
  SDman->AddNewDetector(sensorSD);
}
if(!sampleSD)
{
  sampleSD = new XrayFluoSampleSD ("SampleSD",this);
  SDman->AddNewDetector(sampleSD);
}
if (logicSensor)
{
logicSensor->SetSensitiveDetector(sensorSD);
}
if (logicSample)
{
    logicSample->SetSensitiveDetector(sampleSD); 
} 
// Visualization attributes

  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicSensor->SetVisAttributes(simpleBoxVisAtt);
  logicSample->SetVisAttributes(simpleBoxVisAtt);

  //always return the physical World

PrintApparateParameters();
return physiWorld;
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
	 << SampleSizeYZ/cm
	 << " cm * "
	 << SampleSizeYZ/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << sampleMaterial->GetName() 
	 <<G4endl
	 <<"The sensor is a slice  " << SensorThickness/(1.e-6*m) <<  " micron thick"
	 <<G4endl
	 <<"Size: "<< SensorSizeYZ/cm<< " cm * "<<SensorSizeYZ/cm 
	 <<G4endl
	 <<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSampleMaterial(G4String materialChoice)
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

void XrayFluoDetectorConstruction::SetSensorMaterial(G4String materialChoice)
{
 // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);  
  if (pttoMaterial)
     {sensorMaterial = pttoMaterial;
      logicSensor->SetMaterial(pttoMaterial); 
      PrintApparateParameters();
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSampleThickness(G4double val)
{
  // change Sample thickness 
   SampleThickness = val;
} 

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSensorThickness(G4double val)
{
  // change Sensor thickness 
  SensorThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSampleSizeYZ(G4double val)
{
  // change the transverse size of the sample and recompute the world size
  SampleSizeYZ = val;
} 

 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSensorSizeYZ(G4double val)
{
  // change the transverse size of the sensor and recompute the world size
  SensorSizeYZ = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSensorDistance(G4double val)
{
  // change the distance between the sensor and the sample
  Dist = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSensorAzimuth(G4double val)
{
  // change the angular displacement of the sensor from the X-axes
  Theta = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::SetSensorRotation(G4double val)
{
  // change the angle between the sensitive slice and the YZ-plane
  Phi = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
void XrayFluoDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructApparate());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
