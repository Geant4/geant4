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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4StateManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  ComputeCalorParameters();

  // materials
  DefineMaterials();
  SetAbsorberMaterial("G4_Pb");
  SetGapMaterial("G4_lAr");

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
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // use G4-NIST materials data base
  //
  G4NistManager* man = G4NistManager::Instance();
  fDefaultMaterial = man->FindOrBuildMaterial("G4_Galactic");
  man->FindOrBuildMaterial("G4_Pb");
  man->FindOrBuildMaterial("G4_lAr");

  // print table
  //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  ComputeCalorParameters();

  //
  // World
  //
  fSolidWorld = new G4Box("World",                          //its name
    fWorldSizeX/2,fWorldSizeYZ/2,fWorldSizeYZ/2);           //its size

  fLogicWorld = new G4LogicalVolume(fSolidWorld,            //its solid
                                   fDefaultMaterial,        //its material
                                   "World");                //its name

  fPhysiWorld = new G4PVPlacement(nullptr,                  //no rotation
                                 G4ThreeVector(),           //at (0,0,0)
                                 fLogicWorld,               //its logical volume
                                 "World",                   //its name
                                 nullptr,                   //its mother  volume
                                 false,                     //no boolean operation
                                 0);                        //copy number

  //
  // Calorimeter
  //
  fSolidCalor = nullptr; fLogicCalor = nullptr; fPhysiCalor = nullptr;
  fSolidLayer = nullptr; fLogicLayer = nullptr; fPhysiLayer = nullptr;

  if (fCalorThickness > 0.)  {
    fSolidCalor = new G4Box("Calorimeter",                //its name
      fCalorThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);   //size

    fLogicCalor = new G4LogicalVolume(fSolidCalor,        //its solid
                                      fDefaultMaterial,   //its material
                                      "Calorimeter");     //its name

    fPhysiCalor = new G4PVPlacement(nullptr,              //no rotation
                                   G4ThreeVector(),       //at (0,0,0)
                                   fLogicCalor,           //its logical volume
                                   "Calorimeter",         //its name
                                   fLogicWorld,           //its mother  volume
                                   false,                 //no boolean operation
                                   0);                    //copy number

    //
    // Layer
    //
    fSolidLayer = new G4Box("Layer",                      //its name
      fLayerThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);   //size

    fLogicLayer = new G4LogicalVolume(fSolidLayer,        //its solid
                                     fDefaultMaterial,    //its material
                                     "Layer");            //its name
    if (fNbOfLayers > 1) {
      fPhysiLayer = new G4PVReplica("Layer",              //its name
                                   fLogicLayer,           //its logical volume
                                   fLogicCalor,           //its mother
                                   kXAxis,                //axis of replication
                                   fNbOfLayers,           //number of replica
                                   fLayerThickness);      //witdth of replica
    }
    else {
      fPhysiLayer = new G4PVPlacement(nullptr,            //no rotation
                                   G4ThreeVector(),       //at (0,0,0)
                                   fLogicLayer,           //its logical volume
                                   "Layer",               //its name
                                   fLogicCalor,           //its mother  volume
                                   false,                 //no boolean operation
                                   0);                    //copy number
    }
  }

  //
  // Absorber
  //
  fSolidAbsorber = nullptr; fLogicAbsorber = nullptr; fPhysiAbsorber = nullptr;

  if (fAbsorberThickness > 0.) {
    fSolidAbsorber = new G4Box("Absorber",                 //its name
      fAbsorberThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2); //size

    fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    //its solid
                                         fAbsorberMaterial, //its material
                                         fAbsorberMaterial->GetName()); //name

    fPhysiAbsorber = new G4PVPlacement(nullptr,             //no rotation
                                       G4ThreeVector(-fGapThickness/2,0.,0.),  //its position
                                       fLogicAbsorber,      //its logical volume
                                       fAbsorberMaterial->GetName(), //its name
                                       fLogicLayer,         //its mother
                                       false,               //no boulean operat
                                       0);                  //copy number

  }

  //
  // Gap
  //
  fSolidGap = nullptr; fLogicGap = nullptr; fPhysiGap = nullptr;

  if (fGapThickness > 0.) {
    fSolidGap = new G4Box("Gap",
                          fGapThickness/2,fCalorSizeYZ/2,fCalorSizeYZ/2);

    fLogicGap = new G4LogicalVolume(fSolidGap,
                                    fGapMaterial,
                                    fGapMaterial->GetName());

    fPhysiGap = new G4PVPlacement(nullptr,                  //no rotation
                                  G4ThreeVector(fAbsorberThickness/2,0.,0.),   //its position
                                  fLogicGap,               //its logical volume
                                  fGapMaterial->GetName(), //its name
                                  fLogicLayer,             //its mother
                                  false,                   //no boulean operat
                                  0);                      //copy number
  }

  PrintCalorParameters();

  //
  // Visualization attributes
  //
  fLogicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto  simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  fLogicCalor->SetVisAttributes(simpleBoxVisAtt);

  //
  //always return the physical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << fNbOfLayers << " layers of: [ "
         << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName()
         << " + "
         << fGapThickness/mm << "mm of " << fGapMaterial->GetName() << " ] "
         << "\n------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial(const G4String& materialChoice)
{
  // search the material by its name
  auto material = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (material != nullptr) {
    fAbsorberMaterial = material;
    if (fLogicAbsorber != nullptr) {
      fLogicAbsorber->SetMaterial(fAbsorberMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapMaterial(const G4String& materialChoice)
{
  auto material = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
  if (material != nullptr) {
    fGapMaterial = material;
    if ( fLogicGap != nullptr) {
      fLogicGap->SetMaterial(fGapMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberThickness(G4double value)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = value;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGapThickness(G4double value)
{
  // change Gap thickness and recompute the calorimeter parameters
  fGapThickness = value;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double value)
{
  // change the transverse size and recompute the calorimeter parameters
  fCalorSizeYZ = value;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int value)
{
  fNbOfLayers = value;
  if ( G4StateManager::GetStateManager()->GetCurrentState() != G4State_PreInit ) {
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
