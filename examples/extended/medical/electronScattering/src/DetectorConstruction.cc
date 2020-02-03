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
/// \file medical/electronScattering/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:fMaterial_World(0), fMaterial_Frame(0),
 fMaterial_ExitWindow(0), fMaterial_ScatterFoil(0), fMaterial_MonitorChbr(0),
 fMaterial_Bag(0), fMaterial_Gas(0), fMaterial_Ring(0),
 fPvol_World(0), fPvol_Frame(0), fDetectorMessenger(0)
{              
  // materials  
  DefineMaterials();
  
  // geometry
  GeometryParameters();
 
  // create commands for interactive definition of the calorimeter  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{  
G4double a, z, density;
G4int ncomponents, natoms;
G4double fractionmass;
G4double temperature, pressure;

// define Elements
//
G4Element* H  = new G4Element("Hydrogen", "H",  z= 1, a=   1.0079*g/mole);
G4Element* He = new G4Element("Helium",   "He", z= 2, a=   4.0026*g/mole);
G4Element* Be = new G4Element("Beryllium","Be", z= 4, a=   9.1218*g/mole);
G4Element* C  = new G4Element("Carbon",   "C",  z= 6, a=  12.0107*g/mole);
G4Element* N  = new G4Element("Nitrogen", "N",  z= 7, a=  14.0067*g/mole);
G4Element* O  = new G4Element("Oxygen",   "O",  z= 8, a=  15.9994*g/mole);
G4Element* Al = new G4Element("Aluminium","Al", z=13, a=  26.9815*g/mole);
G4Element* Ar = new G4Element("Argon",    "Ar", z=18, a=  39.9480*g/mole);
G4Element* Ti = new G4Element("Titanium", "Ti", z=22, a=  47.8670*g/mole);
G4Element* Va = new G4Element("Vanadium", "Va", z=23, a=  50.9415*g/mole);
G4Element* Cu = new G4Element("Copper",   "Cu", z=29, a=  63.5460*g/mole);
G4Element* Ta = new G4Element("Tantalum", "Ta", z=73, a= 180.9479*g/mole);
G4Element* Au = new G4Element("Gold",     "Au", z=79, a= 196.9666*g/mole);

// Air
//
G4Material* Air = 
new G4Material("Air", density= 1.205*mg/cm3, ncomponents=4,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
Air->AddElement(C, fractionmass=0.000124);                      
Air->AddElement(N, fractionmass=0.755267);
Air->AddElement(O, fractionmass=0.231782);
Air->AddElement(Ar,fractionmass=0.012827);

// Titanium
//
G4Material* Titanium = 
new G4Material("Titanium", density= 4.42*g/cm3, ncomponents=3);
Titanium->AddElement(Ti, fractionmass=0.90);
Titanium->AddElement(Al, fractionmass=0.06);
Titanium->AddElement(Va, fractionmass=0.04);

// Mylar
//
G4Material* Mylar = 
new G4Material("Mylar", density= 1.40*g/cm3, ncomponents=3);
Mylar->AddElement(H, natoms=4);
Mylar->AddElement(C, natoms=5);
Mylar->AddElement(O, natoms=2);

// Helium
//
G4Material* Helium = 
new G4Material("Helium", density= 0.166*mg/cm3, ncomponents=1,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
Helium->AddElement(He, fractionmass=1.0);                      

// Aluminium
//
G4Material* Aluminium = 
new G4Material("Aluminium", density= 2.7*g/cm3, ncomponents=1);
Aluminium->AddElement(Al, fractionmass=1.0);

// Beryllium
//
G4Material* Beryllium = 
new G4Material("Beryllium", density= 1.85*g/cm3, ncomponents=1);
Beryllium->AddElement(Be, fractionmass=1.0);

//Graphite
//
G4Material* Graphite = 
new G4Material("Graphite", density= 2.18*g/cm3, ncomponents=1);
Graphite->AddElement(C, fractionmass=1.0);

// Copper
//
G4Material* Copper = 
new G4Material("Copper", density= 8.92*g/cm3, ncomponents=1);
Copper->AddElement(Cu, fractionmass=1.0);

// Tantalum
//
G4Material* Tantalum = 
new G4Material("Tantalum", density= 16.65*g/cm3, ncomponents=1);
Tantalum->AddElement(Ta, fractionmass=1.0);

// Gold
//
G4Material* Gold = 
new G4Material("Gold", density= 19.30*g/cm3, ncomponents=1);
Gold->AddElement(Au, fractionmass=1.0);

// example of vacuum
//
density     = universe_mean_density;    //from PhysicalConstants.h
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material* Vacuum =
new G4Material("Galactic", z=1, a=1.01*g/mole,density,
                          kStateGas,temperature,pressure);
                          
//print
//
G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
                            
// assign materials
//
fMaterial_World       = Vacuum;
fMaterial_Frame       = Air;
fMaterial_ExitWindow  = Titanium;
fMaterial_ScatterFoil = fMaterial_Frame;
fMaterial_MonitorChbr = Mylar;
fMaterial_Bag         = Mylar;
fMaterial_Gas         = Helium;
fMaterial_Ring        = Aluminium;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::GeometryParameters()
{
  fZfront_ExitWindow = 0.0*um;  
  fThickness_ExitWindow = 41.2*um;
  
  fZfront_ScatterFoil = 2.65*cm;  
  fThickness_ScatterFoil = 0.0*um;
  
  fZfront_MonitorChbr = 50.*mm;  
  fThickness_MonitorChbr = 112.7*um;
  
  fZfront_Bag = 64.975*mm;  
  fThickness_Bag = 110.0050*cm;
  
  fThickness_Gas = 110.*cm;
  
  fThickness_Ring = 14.*mm;
  fInnerRadius_Ring = 20.*cm;
  
  fZfront_Frame = 2.0*um;    
  fThickness_Frame = 118.2*cm;
      
  fThickness_World = fZfront_Frame + fThickness_Frame;
  fRadius_World = 23.3*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{ 
  // Cleanup old geometry
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
        
  // World
  //
  G4Tubs*
  svol_World = new G4Tubs("World",                        //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_World,            //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                           
  lvol_World = new G4LogicalVolume(svol_World,            //its solid
                                   fMaterial_World,       //its material
                                   "World");              //its name
                                   
  fPvol_World = new G4PVPlacement(0,                      //no rotation
                                  G4ThreeVector(),        //no translation
                                  lvol_World,             //its logical volume
                                  "World",                //its name
                                  0,                      //its mother  volume
                                  false,                  //no boolean operation
                                  0);                     //copy number
                                 
  // Frame
  //
  G4Tubs*
  svol_Frame = new G4Tubs("Frame",                        //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_Frame,            //half-length 
                         0., twopi);                      //theta1, theta2
                         
  G4LogicalVolume*                            
  lvol_Frame = new G4LogicalVolume(svol_Frame,            //its solid
                                   fMaterial_Frame,       //its material
                                   "Frame");              //its name
                                   
  G4double 
  zpos = fZfront_Frame;
                                     
  fPvol_Frame = new G4PVPlacement(0,                      //no rotation
                        G4ThreeVector(0,0,zpos),          //translation
                                 lvol_Frame,              //its logical volume
                                 "Frame",                 //its name
                                 lvol_World,              //its mother  volume
                                 false,                   //no boolean operation
                                 0);                      //copy number

                                                
  // ExitWindow
  //
  G4Tubs*
  svol_ExitWindow = new G4Tubs("ExitWindow",              //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_ExitWindow,       //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                         
  lvol_ExitWindow = new G4LogicalVolume(svol_ExitWindow,  //solid
                                   fMaterial_ExitWindow,  //material
                                   "ExitWindow");         //name
                                   
  zpos = fZfront_ExitWindow + 0.5*fThickness_ExitWindow - 0.5*fThickness_Frame;
  
                   new G4PVPlacement(0,                   //no rotation
                            G4ThreeVector(0,0,zpos),      //translation
                                 lvol_ExitWindow,         //logical volume
                                 "ExitWindow",            //name
                                 lvol_Frame,              //mother volume
                                 false,                   //no boolean operation
                                 0);                      //copy number
                                                
  // Monitor Chamber
  //
  G4Tubs*
  svol_MonitorChbr = new G4Tubs("MonitorChbr",            //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_MonitorChbr,      //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                         
  lvol_MonitorChbr = new G4LogicalVolume(svol_MonitorChbr,//solid
                                   fMaterial_MonitorChbr, //material
                                   "MonitorChbr");        //name
                                   
  zpos = fZfront_MonitorChbr + 0.5*fThickness_MonitorChbr - 0.5*fThickness_Frame;
  
                     new G4PVPlacement(0,                 //no rotation
                            G4ThreeVector(0,0,zpos),      //translation
                                 lvol_MonitorChbr,        //logical volume
                                 "MonitorChbr",           //name
                                 lvol_Frame,              //mother volume
                                 false,                   //no boolean operation
                                 0);                      //copy number
                                 
                                                
  // Bag
  //
  G4Tubs*
  svol_Bag = new G4Tubs("Bag",                            //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_Bag,              //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                         
  lvol_Bag = new G4LogicalVolume(svol_Bag,                //solid
                                 fMaterial_Bag,           //material
                                 "Bag");                  //name
                                   
  zpos = fZfront_Bag + 0.5*fThickness_Bag - 0.5*fThickness_Frame;
  
             new G4PVPlacement(0,                         //no rotation
                            G4ThreeVector(0,0,zpos),      //translation
                                 lvol_Bag,                //logical volume
                                 "Bag",                   //name
                                 lvol_Frame,              //mother volume
                                 false,                   //no boolean operation
                                 0);                      //copy number
                                 
                                                
  // Gas
  //
  G4Tubs*
  svol_Gas = new G4Tubs("Gas",                            //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_Gas,              //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                         
  lvol_Gas = new G4LogicalVolume(svol_Gas,                //solid
                                 fMaterial_Gas,           //material
                                 "Gas");                  //name
  

             new G4PVPlacement(0,                         //no rotation
                            G4ThreeVector(),              //no translation
                                 lvol_Gas,                //logical volume
                                 "Gas",                   //name
                                 lvol_Bag,                //mother volume
                                 false,                   //no boolean operation
                                 0);                      //copy number


  // Rings
  //
  G4Tubs*
  svol_Ring = new G4Tubs("Ring",                          //name
                       fInnerRadius_Ring, fRadius_World,  //r1, r2
                         0.5*fThickness_Ring,             //half-length 
                         0., twopi);                      //theta1, theta2

  G4LogicalVolume*                         
  lvol_Ring = new G4LogicalVolume(svol_Ring,              //solid
                                   fMaterial_Ring,        //material
                                   "Ring");               //name
                                   
  zpos = 0.5*(fThickness_Gas - fThickness_Ring);
  
              new G4PVPlacement(0,                        //no rotation
                            G4ThreeVector(0,0,zpos),      //translation
                                 lvol_Ring,               //logical volume
                                 "Ring",                  //name
                                 lvol_Gas,                //mother volume
                                 false,                   //no boolean operation
                                 1);                      //copy number

              new G4PVPlacement(0,                        //no rotation
                            G4ThreeVector(0,0,-zpos),     //translation
                                 lvol_Ring,               //logical volume
                                 "Ring",                  //name
                                 lvol_Gas,                //mother volume
                                 false,                   //no boolean operation
                                 2);                      //copy number

                                                
  // ScatterFoil (only if it is not Air)
  //
  if ((fMaterial_ScatterFoil != fMaterial_Frame) && (fThickness_ScatterFoil > 0.))
   {
    G4Tubs*
    svol_ScatterFoil = new G4Tubs("ScatterFoil",          //name
                         0*cm, fRadius_World,             //r1, r2
                         0.5*fThickness_ScatterFoil,      //half-length 
                         0., twopi);                      //theta1, theta2

    G4LogicalVolume*                         
    lvol_ScatterFoil = new G4LogicalVolume(svol_ScatterFoil,//solid
                                   fMaterial_ScatterFoil,   //material
                                   "ScatterFoil");          //name
                                   
    zpos = fZfront_ScatterFoil + 0.5*fThickness_ScatterFoil - 0.5*fThickness_Frame;
  
                   new G4PVPlacement(0,                  //no rotation
                            G4ThreeVector(0,0,zpos),     //translation
                                 lvol_ScatterFoil,       //logical volume
                                 "ScatterFoil",          //name
                                 lvol_Frame,             //mother volume
                                 false,                  //no boolean operation
                                 0);                     //copy number
  }                                 
                                        
   PrintGeometry();         
  
  //always return the physical World
  //
  return fPvol_World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintGeometry()
{
   
  // choose printing format
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int prec = G4cout.precision(6);
  
  G4cout << "\n \t \t" << "Material \t" << "Z_front \t" << "Thickness \n";  
     
  G4cout << "\n  ExitWindow \t" << fMaterial_ExitWindow->GetName()
         << "\t" << G4BestUnit(fZfront_ExitWindow,"Length")
         << "\t" << G4BestUnit(fThickness_ExitWindow,"Length");
         
  if (fMaterial_ScatterFoil != fMaterial_Frame) {         
    G4cout << "\n  ScatterFoil \t" << fMaterial_ScatterFoil->GetName() << "\t"
           << "\t" << G4BestUnit(fZfront_ScatterFoil,"Length")
           << "\t" << G4BestUnit(fThickness_ScatterFoil,"Length");
  }
                    
  G4cout << "\n  MonitorChbr \t" << fMaterial_MonitorChbr->GetName() << "\t"
         << "\t" << G4BestUnit(fZfront_MonitorChbr,"Length")
         << "\t" << G4BestUnit(fThickness_MonitorChbr,"Length");

  G4double thickBagWindow = 0.5*(fThickness_Bag - fThickness_Gas);
  G4double zfrontGas = fZfront_Bag + thickBagWindow;
  G4double zfrontBagWindow2 = zfrontGas + fThickness_Gas;
            
  G4cout << "\n  BagWindow1 \t" << fMaterial_Bag->GetName() << "\t"
         << "\t" << G4BestUnit(fZfront_Bag,"Length")
         << "\t" << G4BestUnit(thickBagWindow,"Length");
            
  G4cout << "\n  Gas       \t" << fMaterial_Gas->GetName() << "\t"
         << "\t" << G4BestUnit(zfrontGas,"Length")
         << "\t" << G4BestUnit(fThickness_Gas,"Length");
                     
  G4cout << "\n  BagWindow2 \t" << fMaterial_Bag->GetName() << "\t"
         << "\t" << G4BestUnit(zfrontBagWindow2,"Length")
         << "\t" << G4BestUnit(thickBagWindow,"Length");
                     
  G4cout << "\n  ScoringPlane \t" << fMaterial_Frame->GetName() << "\t"
         << "\t" << G4BestUnit(fThickness_Frame,"Length") << "\n" << G4endl;
         
  // restaure default formats
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialScatter(G4String material)
{
  // search the material by its name
  G4Material* pMaterial = G4Material::GetMaterial(material);

  if (pMaterial) fMaterial_ScatterFoil = pMaterial;                  
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetThicknessScatter(G4double val)
{
  fThickness_ScatterFoil = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#include "G4RunManager.hh" 
 
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

