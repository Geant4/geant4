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
/// \file field/field01/src/F01DetectorConstruction.cc
/// \brief Implementation of the F01DetectorConstruction class
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "F01DetectorConstruction.hh"
#include "F01DetectorMessenger.hh"

#include "F01CalorimeterSD.hh"
#include "F01FieldSetup.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01DetectorConstruction::F01DetectorConstruction()
 : fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
   fSolidAbsorber(0),fLogicAbsorber(0), fPhysiAbsorber(0),
   fEmFieldSetup(0), fDetectorMessenger(0), fCalorimeterSD(0),
   fAbsorberMaterial(0), fAbsorberThickness(0.),fAbsorberRadius(0.),
   fWorldChanged(false), fZAbsorber(0.), fZStartAbs(0.), fZEndAbs(0.),
   fWorldMaterial(0), fWorldSizeR(0.), fWorldSizeZ(0.)
{
  // default parameter values of the calorimeter

  fWorldSizeZ = 44000.*mm;
  fWorldSizeR = 22000.*mm;

  fAbsorberThickness = 1.0*mm;

  fAbsorberRadius   = 20000.*mm;

  fZAbsorber = 21990.0*mm ;

  // create commands for interactive definition of the calorimeter  

  fDetectorMessenger = new F01DetectorMessenger(this);  

  // create materials

  DefineMaterials();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

F01DetectorConstruction::~F01DetectorConstruction()
{ 
  delete fDetectorMessenger;
  if (fEmFieldSetup) delete fEmFieldSetup;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* F01DetectorConstruction::Construct()
{
   G4VPhysicalVolume* calorWorld= ConstructCalorimeter();

   // Construct the field creator - this will register the field it creates
   fEmFieldSetup = new F01FieldSetup(G4ThreeVector( 3.3*tesla, 0.0, 0.0 ) ) ; 
   return calorWorld; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
  G4String name, symbol ;             // a=mass of a mole;
  G4double a, z, density ;            // z=mean number of protons;  
  G4int nel ;
  G4int ncomponents;
  G4double fractionmass, pressure, temperature;

  //
  // define Elements
  //

  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  //
  // define simple materials
  //

  // Mylar

  density = 1.39*g/cm3;
  G4Material* Mylar = new G4Material(name="Mylar", density, nel=3);
  Mylar->AddElement(elO,2);
  Mylar->AddElement(elC,5);
  Mylar->AddElement(elH,4);

  // Polypropelene

  G4Material* CH2 = new G4Material ("Polypropelene" , 0.91*g/cm3, 2);
  CH2->AddElement(elH,2);
  CH2->AddElement(elC,1);

  // Krypton as detector gas, STP

  density = 3.700*mg/cm3 ;
  a = 83.80*g/mole ;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  // Dry air (average composition)

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  density  = 1.2928*mg/cm3 ;       // STP
  density *= 1.0e-8 ;       // pumped vacuum

  temperature = STP_Temperature;
  pressure = 1.0e-8*STP_Pressure;

  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3,
                                   kStateGas,temperature,pressure);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

  // Carbon dioxide, STP

  density = 1.842*mg/cm3;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3 ;      
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2", density, ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2", density, ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 ) ;


  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // default materials of the calorimeter
  
  fAbsorberMaterial = Air ; //  Kr20CO2 ;   // XeCO2CF4  ; 

  fWorldMaterial    = Air ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* F01DetectorConstruction::ConstructCalorimeter()
{
  // complete the Calor parameters definition and Print 

  ComputeCalorParameters();
  PrintCalorParameters();
      
  // Cleanup old geometry

  if (fPhysiWorld)
  {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
  }

  // World
  
  fSolidWorld = new G4Tubs("World",                        // its name
                   0.,fWorldSizeR,fWorldSizeZ/2.,0.,twopi);// its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,           // its solid
                                   fWorldMaterial,         // its material
                                   "World");               // its name
                                   
  fPhysiWorld = new G4PVPlacement(0,                       // no rotation
                                  G4ThreeVector(),         // at (0,0,0)
                                  "World",                 // its name
                                  fLogicWorld,             // its logical volume
                                  0,                       // its mother  volume
                                  false,                   // no boolean operation
                                  0);                      // copy number
  // Absorber

  if (fAbsorberThickness > 0.) 
  { 
      fSolidAbsorber = new G4Tubs("Absorber", 1.0*mm, 
                                  fAbsorberRadius,
                                  fAbsorberThickness/2., 
                                  0.0,twopi); 
                          
      fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    
                                           fAbsorberMaterial, 
                                           "Absorber");     
                                                
      fPhysiAbsorber = new G4PVPlacement(0,                   
                                         G4ThreeVector(0.,0.,fZAbsorber),        
                                         "Absorber",        
                                         fLogicAbsorber,     
                                         fPhysiWorld,       
                                         false,             
                                         0);
  }
                                 
  // Sensitive Detectors: Absorber 
  
  G4SDManager* theSDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new F01CalorimeterSD("CalorSD",this);
    theSDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of " 
         << fWorldSizeZ/mm << "mm of " << fWorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << fWorldSizeR/mm << " mm. " 
         << G4endl;
  G4cout << " The ABSORBER is made of " 
         << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName() ;
  G4cout << ", the transverse size (R) is " << fAbsorberRadius/mm << " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " << fZAbsorber/mm << "  mm." 
         << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* material;
  for (size_t j=0 ; j<theMaterialTable->size() ; j++)
   { material = (*theMaterialTable)[j];     
     if (material->GetName() == materialChoice)
        {
          fAbsorberMaterial = material;
          fLogicAbsorber->SetMaterial(material); 
        }             
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* material;
  for (size_t j=0 ; j<theMaterialTable->size() ; j++)
   { material = (*theMaterialTable)[j];     
     if(material->GetName() == materialChoice)
        {
          fWorldMaterial = material;
          fLogicWorld->SetMaterial(material); 
        }             
   }
}

///////////////////////////////////////////////////////////////////////////
//
//

void F01DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = val;
  ComputeCalorParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
//

void F01DetectorConstruction::SetAbsorberRadius(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fAbsorberRadius = val;
  ComputeCalorParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldChanged = true;
  fWorldSizeZ = val;
  ComputeCalorParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::SetWorldSizeR(G4double val)
{
  fWorldChanged = true;
  fWorldSizeR = val;
  ComputeCalorParameters();
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fZAbsorber = val;
  ComputeCalorParameters();
}  


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void F01DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

