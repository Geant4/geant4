// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5DetectorConstruction.cc,v 1.4 2000-01-20 15:34:39 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em5DetectorConstruction.hh"
#include "Em5DetectorMessenger.hh"

#include "Em5CalorimeterSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5DetectorConstruction::Em5DetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidAbsorber(NULL),logicAbsorber(NULL),physiAbsorber(NULL),
 AbsorberMaterial(NULL),WorldMaterial(NULL),
 magField(NULL),calorimeterSD(NULL),defaultWorld(true)
{
  // default parameter values of the calorimeter
  AbsorberThickness = 1.*cm;
  AbsorberSizeYZ    = 2.*cm;
  XposAbs           = 0.*cm ;
  ComputeCalorParameters();

  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new Em5DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em5DetectorConstruction::~Em5DetectorConstruction()
{ 
  delete detectorMessenger;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Em5DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::DefineMaterials()
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
G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

a = 12.01*g/mole;
G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

a = 14.01*g/mole;
G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 16.00*g/mole;
G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 39.95*g/mole;
G4Element* elAr  = new G4Element(name="Argon"  ,symbol="Ar", z= 18., a);

a = 131.29*g/mole;
G4Element* elXe  = new G4Element(name="Xenon"  ,symbol="Xe", z= 54., a);

//
// define simple materials
//
density = 1.848*g/cm3;
a = 9.01*g/mole;
G4Material* Be = new G4Material(name="Beryllium", z=4., a, density);

density = 2.700*g/cm3;
a = 26.98*g/mole;
G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

density = 2.330*g/cm3;
a = 28.09*g/mole;
G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

density = 1.390*g/cm3;
G4Material* lAr = new G4Material(name="liquidArgon", density, ncomponents=1);
lAr->AddElement(elAr, natoms=1);

density = 7.870*g/cm3;
a = 55.85*g/mole;
G4Material* Fe = new G4Material(name="Iron"   , z=26., a, density);

density = 8.960*g/cm3;
a = 63.55*g/mole;
G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

density = 19.32*g/cm3;
a =196.97*g/mole;
G4Material* Au = new G4Material(name="Gold"   , z=79., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//

density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(elH, natoms=2);
H2O->AddElement(elO, natoms=1);

//
// define a material from elements.   case 2: mixture by fractional mass
//

density = 1.290*mg/cm3;
G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
Air->AddElement(elN, fractionmass=0.7);
Air->AddElement(elO, fractionmass=0.3);

//
// examples of gas
//

density = 5.458*mg/cm3;
G4Material* Xe = new G4Material(name="XenonGas" , density, ncomponents=1,
                   kStateGas,temperature=293.15*kelvin, pressure=1*atmosphere);
Xe->AddElement(elXe, natoms=1);

density = 1.977*mg/cm3;
G4Material* CO2 = new G4Material(name="CarbonicGas", density, ncomponents=2);
CO2->AddElement(elC, natoms=1);
CO2->AddElement(elO, natoms=2);

density = 1.8223*mg/cm3;
G4Material* ArCO2 = new G4Material(name="ArgonCO2", density, ncomponents=2);
ArCO2->AddElement (elAr, fractionmass=0.7844);
ArCO2->AddMaterial(CO2,  fractionmass=0.2156);

//another way to define mixture of gas per volume
density = 1.8223*mg/cm3;
G4Material* NewArCO2 = new G4Material(name="NewArgonCO2", density, ncomponents=3);
NewArCO2->AddElement (elAr, natoms=8);
NewArCO2->AddElement (elC,  natoms=2);
NewArCO2->AddElement (elO,  natoms=4);

//Argon-CH4
density = 1.709*mg/cm3;
G4Material* ArCH4 = new G4Material(name="ArgonCH4", density, ncomponents=3);
ArCH4->AddElement (elAr, natoms=93);
ArCH4->AddElement (elC,  natoms=7);
ArCH4->AddElement (elH,  natoms=28);

//Xenon-methane-propane
density = 4.9196*mg/cm3;
G4Material* XeCH = new G4Material(name="XenonMethanePropane", density, ncomponents=3,
                   kStateGas,temperature=293.15*kelvin, pressure=1*atmosphere);
XeCH->AddElement (elXe, natoms=875);
XeCH->AddElement (elC,  natoms=225);
XeCH->AddElement (elH,  natoms=700);

//
// example of vacuum
//

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                     kStateGas,temperature,pressure);
		     
G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //default materials of the calorimeter
  AbsorberMaterial = Si;
  WorldMaterial    = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VPhysicalVolume* Em5DetectorConstruction::ConstructCalorimeter()
{
  // complete the Calor parameters definition and Print 
  ComputeCalorParameters();

  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeX/2,WorldSizeYZ/2,WorldSizeYZ/2);   //its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   WorldMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  
  //                               
  // Absorber
  // 
  solidAbsorber = new G4Box("Absorber",	
                      AbsorberThickness/2,AbsorberSizeYZ/2,AbsorberSizeYZ/2); 
                          
  logicAbsorber = new G4LogicalVolume(solidAbsorber,    //its solid
    	                  AbsorberMaterial,             //its material
   	                 "Absorber");                   //its name
      			                  
  physiAbsorber = new G4PVPlacement(0,		   //no rotation
      		  G4ThreeVector(XposAbs,0.,0.),    //its position
                                "Absorber",        //its name
                                logicAbsorber,     //its logical volume
                                physiWorld,        //its mother
                                false,             //no boulean operat
                                0);                //copy number
                                        
  
  //                               
  // Sensitive Detectors: Absorber 
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!calorimeterSD)
  {
    calorimeterSD = new Em5CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( calorimeterSD );
  }
  if (logicAbsorber)
      logicAbsorber->SetSensitiveDetector(calorimeterSD);
      
  //                                        
  // Visualization attributes
  //
  G4VisAttributes* VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(VisAtt);
  
  //
  //always return the physical World
  //
  PrintCalorParameters();  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of " 
         << G4BestUnit(WorldSizeX,"Length") << " of " << WorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is " 
         << G4BestUnit(WorldSizeYZ,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(AbsorberThickness,"Length")<< " of " << AbsorberMaterial->GetName();
  G4cout << ". The transverse size (YZ) is " 
         << G4BestUnit(AbsorberSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the absorber " << G4BestUnit(XposAbs,"Length");
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {AbsorberMaterial = pttoMaterial;
      logicAbsorber->SetMaterial(pttoMaterial); 
      PrintCalorParameters();
     }                  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {WorldMaterial = pttoMaterial;
      logicWorld->SetMaterial(pttoMaterial); 
      PrintCalorParameters();     
     }
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  AbsorberSizeYZ = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetWorldSizeX(G4double val)
{
  WorldSizeX = val;
  defaultWorld = false;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  WorldSizeYZ = val;
  defaultWorld = false;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetAbsorberXpos(G4double val)
{
  XposAbs  = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em5DetectorConstruction::SetMagField(G4double fieldValue)
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
  
void Em5DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

