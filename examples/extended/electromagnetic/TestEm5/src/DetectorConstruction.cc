//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: DetectorConstruction.cc,v 1.3 2004/06/21 10:57:13 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:AbsorberMaterial(0),WorldMaterial(0),defaultWorld(true),
 solidWorld(0),logicWorld(0),physiWorld(0),
 solidAbsorber(0),logicAbsorber(0),physiAbsorber(0),
 magField(0)
{
  // default parameter values of the calorimeter
  AbsorberThickness = 1.*cm;
  AbsorberSizeYZ    = 2.*cm;
  XposAbs           = 0.*cm;
  ComputeCalorParameters();
  
  // materials  
  DefineMaterials();
  SetWorldMaterial   ("Galactic");
  SetAbsorberMaterial("Silicon");
 
  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String symbol;             //a=mass of a mole;
G4double a, z, density;      //z=mean number of protons;  

G4int ncomponents, natoms;
G4double fractionmass;
G4double temperature, pressure;

//
// define Elements
//

G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.01*g/mole);
G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.01*g/mole);
G4Element* N  = new G4Element("Nitrogen",symbol="N",  z= 7, a=  14.01*g/mole);
G4Element* O  = new G4Element("Oxygen",  symbol="O",  z= 8, a=  16.00*g/mole);
G4Element* Ar = new G4Element("Argon",   symbol="Ar", z=18, a=  39.95*g/mole);
G4Element* Xe = new G4Element("Xenon",   symbol="Xe", z=54, a= 131.29*g/mole);

//
// define simple materials
//

new G4Material("Beryllium", z= 4, a= 9.01*g/mole, density= 1.848*g/cm3);
new G4Material("Aluminium", z=13, a=26.98*g/mole, density= 2.700*g/cm3);
new G4Material("Silicon"  , z=14, a=28.09*g/mole, density= 2.330*g/cm3);

G4Material* lAr = 
new G4Material("liquidArgon", density= 1.390*g/cm3, ncomponents=1);
lAr->AddElement(Ar, natoms=1);

new G4Material("Iron",     z=26, a= 55.85*g/mole, density= 7.870*g/cm3);
new G4Material("Copper",   z=29, a= 63.55*g/mole, density= 8.960*g/cm3);
new G4Material("Tungsten", z=74, a=183.85*g/mole, density= 19.30*g/cm3);
new G4Material("Gold",     z=79, a=196.97*g/mole, density= 19.32*g/cm3);
new G4Material("Lead",     z=82, a=207.19*g/mole, density= 11.35*g/cm3);

//
// define a material from elements.   case 1: chemical molecule
//

G4Material* H2O = new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);
H2O->GetIonisation()->SetMeanExcitationEnergy(75*eV);

//
// define a material from elements.   case 2: mixture by fractional mass
//

G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// examples of gas
//

new G4Material("XenonGas", z=54, a=131.29*g/mole, density= 5.458*mg/cm3,
                           kStateGas,293.15*kelvin, 1*atmosphere);

G4Material* CO2 =
new G4Material("CarbonicGas", density= 1.977*mg/cm3, ncomponents=2);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);

G4Material* ArCO2 =
new G4Material("ArgonCO2",   density= 1.8223*mg/cm3, ncomponents=2);
ArCO2->AddElement (Ar,  fractionmass=0.7844);
ArCO2->AddMaterial(CO2, fractionmass=0.2156);

//another way to define mixture of gas per volume
G4Material* NewArCO2 =
new G4Material("NewArgonCO2", density= 1.8223*mg/cm3, ncomponents=3);
NewArCO2->AddElement (Ar, natoms=8);
NewArCO2->AddElement (C,  natoms=2);
NewArCO2->AddElement (O,  natoms=4);

G4Material* ArCH4 = 
new G4Material("ArgonCH4",    density= 1.709*mg/cm3,  ncomponents=3);
ArCH4->AddElement (Ar, natoms=93);
ArCH4->AddElement (C,  natoms=7);
ArCH4->AddElement (H,  natoms=28);

G4Material* XeCH = 
new G4Material("XenonMethanePropane", density= 4.9196*mg/cm3, ncomponents=3,
                                      kStateGas, 293.15*kelvin, 1*atmosphere);
XeCH->AddElement (Xe, natoms=875);
XeCH->AddElement (C,  natoms=225);
XeCH->AddElement (H,  natoms=700);

//
// example of vacuum
//

density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                                              kStateGas,temperature,pressure);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  xstartAbs = XposAbs-0.5*AbsorberThickness; 
  xendAbs   = XposAbs+0.5*AbsorberThickness;
     
  if (defaultWorld) {
     WorldSizeX = 1.5*AbsorberThickness; WorldSizeYZ= 1.2*AbsorberSizeYZ;
  } 	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{ 
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // complete the Calor parameters definition 
  ComputeCalorParameters();
        
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeX/2,WorldSizeYZ/2,WorldSizeYZ/2);   //its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   WorldMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
                                 
  // Absorber
  // 
  solidAbsorber = new G4Box("Absorber",	
                      AbsorberThickness/2,AbsorberSizeYZ/2,AbsorberSizeYZ/2); 
                          
  logicAbsorber = new G4LogicalVolume(solidAbsorber,    //its solid
    	                  	      AbsorberMaterial, //its material
   	                  	     "Absorber");       //its name
      			                  
  physiAbsorber = new G4PVPlacement(0,		   //no rotation
      		  G4ThreeVector(XposAbs,0.,0.),    //its position
                                logicAbsorber,     //its logical volume
				"Absorber",         //its name
                                logicWorld,        //its mother
                                false,             //no boulean operat
                                0);                //copy number
                                        
   PrintCalorParameters();         
  
  //always return the physical World
  //
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n" << WorldMaterial    << G4endl;
  G4cout << "\n" << AbsorberMaterial << G4endl;
    
  G4cout << "\n The  WORLD   is made of "  << G4BestUnit(WorldSizeX,"Length")
         << " of " << WorldMaterial->GetName();
  G4cout << ". The transverse size (YZ) of the world is " 
         << G4BestUnit(WorldSizeYZ,"Length") << G4endl;
  G4cout << " The ABSORBER is made of " 
         <<G4BestUnit(AbsorberThickness,"Length")
	 << " of " << AbsorberMaterial->GetName();
  G4cout << ". The transverse size (YZ) is " 
         << G4BestUnit(AbsorberSizeYZ,"Length") << G4endl;
  G4cout << " X position of the middle of the absorber "
         << G4BestUnit(XposAbs,"Length");
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) AbsorberMaterial = pttoMaterial;                  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial) WorldMaterial = pttoMaterial;
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  AbsorberThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberSizeYZ(G4double val)
{
  AbsorberSizeYZ = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeX(G4double val)
{
  WorldSizeX = val;
  defaultWorld = false;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeYZ(G4double val)
{
  WorldSizeYZ = val;
  defaultWorld = false;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberXpos(G4double val)
{
  XposAbs  = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void DetectorConstruction::SetMagField(G4double fieldValue)
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

