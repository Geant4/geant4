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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fAbsorMaterial(0), fLAbsor(0), fContainMaterial(0), fLContain(0),
 fWorldMaterial(0), fPWorld(0), fDetectorMessenger(0)
{
  fAbsorRadius = 15*mm;
  fAbsorLength = 60*mm;
  fContainThickness = 2.4*mm;
  DefineMaterials();
  SetAbsorMaterial  ("BeO");
  SetContainMaterial("Stainless-Steel");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4int ncomponents, natoms;
  
  G4Element* Be = new G4Element("Beryllium","Be" ,  4.,  9.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N"  ,  7., 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O"  ,  8., 16.00*g/mole);
  G4Element* Cr = new G4Element("Chromium" ,"Cr" , 24., 51.99*g/mole);
  G4Element* Fe = new G4Element("Iron"     ,"Fe" , 26., 55.84*g/mole);
  G4Element* Ni = new G4Element("Nickel"   ,"Ni" , 28., 58.69*g/mole);
  
  G4Material* BeO = 
  new G4Material("BeO", 3.05*g/cm3, ncomponents=2);
  BeO->AddElement(Be, natoms=1);
  BeO->AddElement( O, natoms=1);
  
  G4Material* inox = 
  new G4Material("Stainless-Steel", 8*g/cm3, ncomponents=3);
  inox->AddElement(Fe, 74*perCent);
  inox->AddElement(Cr, 18*perCent);
  inox->AddElement(Ni,  8*perCent);

  G4Material* Air = 
  new G4Material("Air", 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  fWorldMaterial = Air;
  
  ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  // compute dimensions
  G4double ContainRadius = fAbsorRadius + fContainThickness;
  G4double ContainLength = fAbsorLength + 2*fContainThickness;
  
  G4double WorldSizeXY = 2.4*ContainRadius;
  G4double WorldSizeZ  = 1.2*ContainLength;
  
  // World
  //
  G4Box*
  sWorld = new G4Box("World",                                    //name
              0.5*WorldSizeXY,0.5*WorldSizeXY,0.5*WorldSizeZ);   //dimensions
                   
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMaterial,            //material
                             "World");                  //name

  fPWorld = new G4PVPlacement(0,                        //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number

  // Container
  //
  G4Tubs* 
  sContain = new G4Tubs("Container",                            //name
             0., ContainRadius, 0.5*ContainLength, 0., twopi);  //dimensions

  fLContain = new G4LogicalVolume(sContain,            //shape
                       fContainMaterial,               //material
                       fContainMaterial->GetName());   //name

           new G4PVPlacement(0,                        //no rotation
                       G4ThreeVector(),                //at (0,0,0)
                       fLContain,                      //logical volume
                       fContainMaterial->GetName(),    //name
                       lWorld,                         //mother  volume
                       false,                          //no boolean operation
                       0);                             //copy number

  // Absorber
  //
  G4Tubs* 
  sAbsor = new G4Tubs("Absorber",                                //name
               0., fAbsorRadius, 0.5*fAbsorLength, 0., twopi);    //dimensions

  fLAbsor = new G4LogicalVolume(sAbsor,                //shape
                       fAbsorMaterial,                 //material
                       fAbsorMaterial->GetName());     //name

           new G4PVPlacement(0,                        //no rotation
                       G4ThreeVector(),                //at (0,0,0)
                       fLAbsor,                        //logical volume
                       fAbsorMaterial->GetName(),      //name
                       fLContain,                      //mother  volume
                       false,                          //no boolean operation
                       0);                             //copy number

  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
 G4cout << "\n The Absorber  is a cylinder of " << fAbsorMaterial->GetName()
        << "  radius = " << G4BestUnit(fAbsorRadius,"Length")
        << "  length = " << G4BestUnit(fAbsorLength,"Length") << G4endl;
 G4cout << " The Container is a cylinder of " << fContainMaterial->GetName()
        << "  thickness = " << G4BestUnit(fContainThickness,"Length") << G4endl;
 
 G4cout << "\n" << fAbsorMaterial << G4endl;
 G4cout << "\n" << fContainMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fAbsorMaterial = pttoMaterial;
    if(fLAbsor) { fLAbsor->SetMaterial(fAbsorMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetAbsorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fContainMaterial = pttoMaterial;
    if(fLContain) { fLContain->SetMaterial(fContainMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetContainMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorRadius(G4double value)
{
  fAbsorRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorLength(G4double value)
{
  fAbsorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainThickness(G4double value)
{
  fContainThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

