// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3DetectorConstruction.cc,v 1.1 1999-10-11 16:55:51 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em3DetectorConstruction.hh"
#include "Em3DetectorMessenger.hh"

#include "Em3CalorimeterSD.hh"

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

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3DetectorConstruction::Em3DetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidCalor(NULL),logicCalor(NULL),physiCalor(NULL),
 solidLayer(NULL),logicLayer(NULL),physiLayer(NULL),
 defaultMaterial(NULL),
 magField(NULL),calorimeterSD(NULL)
{
  for (G4int i=0; i<MaxAbsor; i++)
     {
      AbsorMaterial[i]=NULL; AbsorThickness[i]=0.;
      solidAbsor[i]=NULL;logicAbsor[i]=NULL;physiAbsor[i]=NULL;
     }
  //    
  // default parameter values of the calorimeter
  NbOfAbsor = 2;
  AbsorThickness[0] = 2.3*mm;
  AbsorThickness[1] = 5.7*mm;
  NbOfLayers        = 50;
  CalorSizeYZ       = 40.*cm;
  ComputeCalorParameters();

  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new Em3DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em3DetectorConstruction::~Em3DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Em3DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol;             //a=mass of a mole;
G4double a, z, density;            //z=mean number of protons;  
G4int iz, n;                       //iz=number of protons  in an isotope; 
                                   // n=number of nucleons in an isotope;

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

a = 14.01*g/mole;
G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 16.00*g/mole;
G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 28.09*g/mole;
G4Element* Si = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element(name="enriched Uranium", symbol="U", ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//

density = 2.700*g/cm3;
a = 26.98*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

density = 1.390*g/cm3;
a = 39.95*g/mole;
G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

density = 7.700*g/cm3;
a = 55.85*g/mole;
G4Material* Fe = new G4Material(name="Iron"     , z=26., a, density);

density = 8.960*g/cm3;
a = 63.55*g/mole;
G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

density = 19.30*g/cm3;
a = 183.85*g/mole;
G4Material* W = new G4Material(name="Tungsten"  , z=74., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

density = 18.95*g/cm3;
a = 238.03*g/mole;
G4Material* Ur = new G4Material(name="Uranium"  , z=92., a, density);

//
// define a material from elements.   case 1: chemical molecule
//
 
density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);

density = 1.032*g/cm3;
G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
Sci->AddElement(C, natoms=9);
Sci->AddElement(H, natoms=10);

density = 2.330*g/cm3;
G4Material* Sili = new G4Material(name="Silicium", density, ncomponents=1);
Sili->AddElement(Si, natoms=1);

density = 2.200*g/cm3;
G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

density = 1.700*g/cm3;
G4Material* G10 = new G4Material(name="NemaG10", density, ncomponents=4);
G10->AddElement(Si, natoms=1);
G10->AddElement(O , natoms=2);
G10->AddElement(C , natoms=3);
G10->AddElement(H , natoms=3);


//
// define a material from elements.   case 2: mixture by fractional mass
//

density = 1.290*mg/cm3;
G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

density = 0.200*g/cm3;
G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
Aerog->AddMaterial(SiO2, fractionmass=62.5*perCent);
Aerog->AddMaterial(H2O , fractionmass=37.4*perCent);
Aerog->AddElement (C   , fractionmass= 0.1*perCent);

//
// examples of gas in non STP conditions
//

density     = 27.*mg/cm3;
pressure    = 50.*atmosphere;
temperature = 325.*kelvin;
G4Material* CO2 = new G4Material(name="CarbonicGas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
density     = 0.3*mg/cm3;
pressure    = 2.*atmosphere;
temperature = 500.*kelvin;
G4Material* steam = new G4Material(name="WaterSteam", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material* vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

density     = 1.e-5*g/cm3;
pressure    = 2.e-2*bar;
temperature = STP_Temperature;         //from PhysicalConstants.h
G4Material* beam = new G4Material(name="Beam", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
beam->AddMaterial(Air, fractionmass=1.);

G4cout << *(G4Material::GetMaterialTable()) << endl;

  //default materials of the calorimeter
  AbsorMaterial[0] = Pb;
  AbsorMaterial[1] = lAr;
  defaultMaterial  = vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VPhysicalVolume* Em3DetectorConstruction::ConstructCalorimeter()
{
  // complete the Calor parameters definition 
  ComputeCalorParameters();
   
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
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
  
  //                               
  // Calorimeter
  //  
  solidCalor=NULL; logicCalor=NULL; physiCalor=NULL;
  
  solidCalor = new G4Box("Calorimeter",				     //its name
    		       CalorThickness/2,CalorSizeYZ/2,CalorSizeYZ/2);//size
    			     
  logicCalor = new G4LogicalVolume(solidCalor,		//its solid
      				   defaultMaterial,	//its material
      				   "Calorimeter");	//its name
    				       
  physiCalor = new G4PVPlacement(0,			//no rotation
                                 G4ThreeVector(),	//at (0,0,0)
                                 "Calorimeter",		//its name
                                 logicCalor,		//its logical volume
                                 physiWorld,		//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                                 
  // Layers
  //
  solidLayer=NULL; logicLayer=NULL; physiLayer=NULL;
    
  solidLayer = new G4Box("Layer",		                      //its name
                       LayerThickness/2,CalorSizeYZ/2,CalorSizeYZ/2); //size
                       
  logicLayer = new G4LogicalVolume(solidLayer,		//its solid
                                   defaultMaterial,	//its material
                                   "Layer");		//its name
  if (NbOfLayers > 1)                                      
    physiLayer = new G4PVReplica("Layer",		//its name
      				 logicLayer,		//its logical volume
      				 physiCalor,		//its mother
                                 kXAxis,		//axis of replication
                                 NbOfLayers,		//number of replica
                                 LayerThickness);	//witdth of replica
  else
    physiLayer = new G4PVPlacement(0,			//no rotation
                                   G4ThreeVector(),	//at (0,0,0)
                                   "Layer",		//its name
                                   logicLayer,		//its logical volume
                                   physiCalor,		//its mother  volume
                                   false,		//no boolean operation
                                   0);			//copy number     				       

  //                               
  // Absorbers
  //
  for (G4int j=0; j<MaxAbsor; j++)  
     { solidAbsor[j]=NULL; logicAbsor[j]=NULL; physiAbsor[j]=NULL;}  

  G4double xfront = -0.5*LayerThickness;  
  for (G4int k=0; k<NbOfAbsor; k++)
     { solidAbsor[k] = new G4Box("Absorber",		//its name
                          AbsorThickness[k]/2,CalorSizeYZ/2,CalorSizeYZ/2); 
                          
      logicAbsor[k] = new G4LogicalVolume(solidAbsor[k],    //its solid
      			                  AbsorMaterial[k], //its material
      			                  "Absorber");      //its name

      G4double xcenter = xfront+0.5*AbsorThickness[k];
      xfront += AbsorThickness[k];       			                  
      physiAbsor[k] = new G4PVPlacement(0,		   //no rotation
      		    G4ThreeVector(xcenter,0.,0.),          //its position
                                        "Absorber",        //its name
                                        logicAbsor[k],     //its logical volume
                                        physiLayer,        //its mother
                                        false,             //no boulean operat
                                        k);                //copy number
                                        
     }
   
  //                               
  // Sensitive Detectors: all Absorbers
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!calorimeterSD)
  {
    calorimeterSD = new Em3CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( calorimeterSD );
  }
  for (G4int l=0; l<NbOfAbsor; l++)  
      logicAbsor[l]->SetSensitiveDetector(calorimeterSD);
  
  //                                        
  // Visualization attributes
  //
  ////logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicCalor->SetVisAttributes(simpleBoxVisAtt);
  
  //
  //always return the physical World
  //
  PrintCalorParameters();
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << NbOfLayers << " layers of:";
  for (G4int i=0; i<NbOfAbsor; i++)
     { 
       G4cout << "\n \t" << setw(12) << AbsorMaterial[i]->GetName() <<": "
              << setw(6) << G4BestUnit(AbsorThickness[i],"Length");
     }
  G4cout << "\n-------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > MaxAbsor)
    { G4cout << "\n ---> warning from SetNbOfAbsor: " 
             << ival << " must be at least 1 and and most " << MaxAbsor 
	     << ". Command refused" << endl;
      return;
    }  
  NbOfAbsor = ival;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetAbsorMaterial(G4int ival,G4String material)
{
  // search the material by its name
  //
  if (ival >= NbOfAbsor)
    { G4cout << "\n --->warning from SetAbsorMaterial: absor number " 
             << ival << " out of range. Command refused" << endl;
      return;
    }
         
  G4Material* pttoMaterial = G4Material::GetMaterial(material);     
  if (pttoMaterial) AbsorMaterial[ival] = pttoMaterial;             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetAbsorThickness(G4int ival,G4double val)
{
  // change Absorber thickness
  //
  if (ival >= NbOfAbsor)
    { G4cout << "\n --->warning from SetAbsorThickness: absor number " 
             << ival << " out of range. Command refused" << endl;
      return;
    }
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetAbsorThickness: thickness " 
             << val  << " out of range. Command refused" << endl;
      return;
    }     	      
  AbsorThickness[ival] = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN)
    { G4cout << "\n --->warning from SetCalorSizeYZ: thickness " 
             << val  << " out of range. Command refused" << endl;
      return;
    }       
  CalorSizeYZ = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 1)
    { G4cout << "\n --->warning from SetNbOfLayers: " 
             << ival << " must be at least 1. Command refused" << endl;
      return;
    }  
  NbOfLayers = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em3DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  //
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
  
void Em3DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
