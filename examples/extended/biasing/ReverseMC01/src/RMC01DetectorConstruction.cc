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
/// \file biasing/ReverseMC01/src/RMC01DetectorConstruction.cc
/// \brief Implementation of the RMC01DetectorConstruction class
//
//
//////////////////////////////////////////////////////////////
//      Class Name:        RMC01DetectorConstruction
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RMC01DetectorConstruction.hh"
#include "RMC01DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "RMC01SD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorConstruction::RMC01DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fDetectorMessenger(0),
   fShield_Thickness(5.*mm),
   fSensitive_cylinder_H (1.*mm),
   fSensitive_cylinder_Rout (1.*mm)
{ 
   fDetectorMessenger = new RMC01DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorConstruction::~RMC01DetectorConstruction()
{ delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RMC01DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructSimpleGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorConstruction::DefineMaterials()
{ 
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  
  G4double fractionmass;
  G4int ncomponents;
  
  //
  // define Elements
  //
  
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
   
  //
  // define simple materials
  //

  new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.700*g/cm3);
  new G4Material("Silicon", z=14., a=28.09*g/mole, density=2.33*g/cm3);
  new G4Material("Tantalum", z=73., a=180.9479*g/mole, density=16.654*g/cm3);

  //
  // define air   
  //

  G4Material* Air = new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  //
  //Example of Vacuum
  //

   new G4Material("Vacuum", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 3.e-18*pascal, 2.73*kelvin);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RMC01DetectorConstruction::ConstructSimpleGeometry()
{

  // Clean old geometry, if any
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // World
  //-----------
  
  G4Box* solidWorld = new G4Box("World",15.*cm, 15.*cm, 15.*cm);        
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,                
                                   G4Material::GetMaterial("Vacuum"),        
                                   "World");                
                                   
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                 logicWorld,                //its logical volume
                 "World",                //its name
                 0,                        //its mother  volume
                 false,                        //no boolean operation
                 0);
                                 
  //Shielding Aluminum Sphere
  //-------------------
  
  G4double radiusShieldingSphere =10.*cm;
    
  G4Orb* solidShieldingSphere=new G4Orb("Shielding", radiusShieldingSphere);
  G4LogicalVolume* logicShieldingSphere=
                       new G4LogicalVolume(solidShieldingSphere,
                                          G4Material::GetMaterial("Aluminum"),
                                          "Shielding");        //its name;
    
  new G4PVPlacement(0,                        //no rotation
                     G4ThreeVector(),        //at (0,0,0)
                     logicShieldingSphere,        //its logical volume
                     "Shielding",        //its name
                      logicWorld,        //its mother  volume
                     false,                //no boolean operation
                     0);
                                     
  //Bulk Sphere
  //-------------------
   
  G4Orb* solidBulkSphere=new G4Orb("Bulk",
                    radiusShieldingSphere-fShield_Thickness);
  G4LogicalVolume* logicBulkSphere=new G4LogicalVolume(
        solidBulkSphere,//its solid
                  G4Material::GetMaterial("Air"),//its material
                                   "Bulk");        //its name;
    
  new G4PVPlacement(0,                        //no rotation
                      G4ThreeVector(),        //at (0,0,0)
                      logicBulkSphere,        //its logical volume
                      "Bulk",        //its name
                      logicShieldingSphere,        //its mother  volume
                      false,                //no boolean operation
                      0);
   
  //Detecting cylinder
  //-------------------
    
  G4Tubs* solidDetecting=new G4Tubs("SensitiveVolume",
                         0.,fSensitive_cylinder_Rout,fSensitive_cylinder_H/2.,
                                    0.,twopi);
    
  G4LogicalVolume* logicDetectingCylinder=new G4LogicalVolume(solidDetecting,
                                           G4Material::GetMaterial("Silicon"),
                                           "SensitiveVolume");
    
  new G4PVPlacement(0,                        //no rotation
                      G4ThreeVector(0.,0.,0.),        //at (0,0,0)
                      logicDetectingCylinder,        //its logical volume
                      "SensitiveVolume",        //its name
                      logicBulkSphere,        //its mother  volume
                      false,                //no boolean operation
                      0);
                                     
    
  RMC01SD* theSensitiveDetector  = new RMC01SD("/SensitiveCylinder");
    
  G4SDManager::GetSDMpointer()->AddNewDetector(theSensitiveDetector);
  logicDetectingCylinder->SetSensitiveDetector(theSensitiveDetector);
                                     
  //Tantalum Plates on the top and beside
  //-------------------------------------
  G4Box* solidPlate=new G4Box("TantalumPlate",4.*cm,4.*cm,0.25*mm);
  G4LogicalVolume* logicPlate=new G4LogicalVolume(solidPlate,        //its solid
                            G4Material::GetMaterial("Tantalum"),//its material
                                          "TantalumPlate");        //its name;
   
   
  new G4PVPlacement(0,                        //no rotation
                      G4ThreeVector(0.,0.,6.*cm),        //at (0,0,0)
                      logicPlate,        //its logical volume
                      "TantalumPlate1",        //its name
                      logicBulkSphere,        //its mother  volume
                      false,                //no boolean operation
                      0);
   
  new G4PVPlacement(0,                        //no rotation
                   G4ThreeVector(0.,0.,-6.*cm),        //at (0,0,0)
                   logicPlate,        //its logical volume
                   "TantalumPlate2",        //its name
                   logicBulkSphere,        //its mother  volume
                                     false,          //no boolean operation
                    0);
   
  return physiWorld;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorConstruction::SetSensitiveVolumeRadius(G4double r)
{  fSensitive_cylinder_Rout=r;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorConstruction::SetSensitiveVolumeHeight(G4double h)
{  fSensitive_cylinder_H=h;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorConstruction::SetShieldingThickness(G4double d)
{ fShield_Thickness=d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
