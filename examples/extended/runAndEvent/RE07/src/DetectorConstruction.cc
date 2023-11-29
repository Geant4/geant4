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
/// \file src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorMessenger.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : fWorldMaterial(nullptr),
    fLogicWorld(nullptr),
    fPhysiWorld(nullptr),
    fLogicLayerFront(nullptr),
    fLogicLayerBack(nullptr)
{
  for (G4int i = 0; i < kMaxAbsor; ++i) {
    fAbsorMaterial[i] = nullptr;
    fAbsorThickness[i] = 0.0;
    fLogicAbsorFront[i] = nullptr;
    fLogicAbsorBack[i] = nullptr;
  }

  // default parameter values of the calorimeter
  fNbOfAbsor = 2;
  fAbsorThickness[1] = 2.3 * mm;
  fAbsorThickness[2] = 5.7 * mm;
  fNbOfLayers = 50;
  fCalorSizeYZ = 40. * cm;
  ComputeCalorParameters();

  // materials
  SetWorldMaterial("G4_Galactic");
  SetAbsorMaterial(1, "G4_Pb");
  SetAbsorMaterial(2, "G4_lAr");

  // create commands for interactive definition of the calorimeter
  fDetectorMessenger.reset(new DetectorMessenger(this));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ComputeCalorParameters()
{
  // Compute derived parameters of the calorimeter
  fLayerThickness = 0.;
  for (G4int iAbs = 1; iAbs <= fNbOfAbsor; iAbs++) {
    fLayerThickness += fAbsorThickness[iAbs];
  }
  fCalorThickness = fNbOfLayers * fLayerThickness;
  fWorldSizeX = 1.2 * fCalorThickness;
  fWorldSizeYZ = 1.2 * fCalorSizeYZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (fPhysiWorld) {
    return fPhysiWorld;
  }
  // complete the Calor parameters definition
  ComputeCalorParameters();

  //
  // World
  //
  auto* solidWorld = new G4Box("World",  // its name
    fWorldSizeX / 2, fWorldSizeYZ / 2,
    fWorldSizeYZ / 2);  // its size

  fLogicWorld = new G4LogicalVolume(solidWorld,  // its solid
    fWorldMaterial,  // its material
    "World");  // its name

  fPhysiWorld = new G4PVPlacement(0,  // no rotation
    G4ThreeVector(),  // at (0,0,0)
    fLogicWorld,  // its fLogical volume
    "World",  // its name
    0,  // its mother  volume
    false,  // no boolean operation
    0);  // copy number
  //
  // Calorimeter
  //

  auto* solidCalor =
    new G4Box("Calorimeter", fCalorThickness / 2, fCalorSizeYZ / 2, fCalorSizeYZ / 2);

  auto* logicCalor = new G4LogicalVolume(solidCalor, fWorldMaterial, "Calorimeter");

  new G4PVPlacement(0,  // no rotation
    G4ThreeVector(),  // at (0,0,0)
    logicCalor,  // its fLogical volume
    "Calorimeter",  // its name
    fLogicWorld,  // its mother  volume
    false,  // no boolean operation
    0);  // copy number

  //
  // Layers
  //

  auto* solidLayer = new G4Box("Layer", fLayerThickness / 2, fCalorSizeYZ / 2, fCalorSizeYZ / 2);

  fLogicLayerFront = new G4LogicalVolume(solidLayer, fWorldMaterial, "Layer-front");
  fLogicLayerBack = new G4LogicalVolume(solidLayer, fWorldMaterial, "Layer-back");
  G4double xfront = -0.5 * fCalorThickness;
  for (G4int l = 0; l < fNbOfLayers; ++l) {
    G4double xcenter = xfront + 0.5 * fLayerThickness;
    xfront += fLayerThickness;
    G4LogicalVolume* logicLayer = fLogicLayerFront;
    if (xcenter > 0) {
      logicLayer = fLogicLayerBack;
    }

    new G4PVPlacement(0, G4ThreeVector(xcenter, 0, 0), logicLayer, "Layer", logicCalor, false, l);
  }

  //
  // Regions
  //

  auto* regionFront = new G4Region("Front");
  regionFront->SetProductionCuts(
    G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
  regionFront->AddRootLogicalVolume(fLogicLayerFront);
  auto* regionBack = new G4Region("Back");
  regionBack->SetProductionCuts(
    G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
  regionBack->AddRootLogicalVolume(fLogicLayerBack);

  //
  // Absorbers
  //

  xfront = -0.5 * fLayerThickness;
  for (G4int k = 1; k <= fNbOfAbsor; ++k) {
    auto* solidAbsor = new G4Box("Absorber",  // its name
      fAbsorThickness[k] / 2, fCalorSizeYZ / 2, fCalorSizeYZ / 2);

    fLogicAbsorFront[k] = new G4LogicalVolume(solidAbsor,  // its solid
      fAbsorMaterial[k],  // its material
      fAbsorMaterial[k]->GetName());
    fLogicAbsorBack[k] = new G4LogicalVolume(solidAbsor,  // its solid
      fAbsorMaterial[k],  // its material
      fAbsorMaterial[k]->GetName());

    G4double xcenter = xfront + 0.5 * fAbsorThickness[k];
    xfront += fAbsorThickness[k];
    new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicAbsorFront[k],
      fAbsorMaterial[k]->GetName(), fLogicLayerFront, false,
      k);  // copy number
    new G4PVPlacement(0, G4ThreeVector(xcenter, 0., 0.), fLogicAbsorBack[k],
      fAbsorMaterial[k]->GetName(), fLogicLayerBack, false,
      k);  // copy number
  }

  PrintCalorParameters();

  // always return the fPhysical World
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n-------------------------------------------------------------"
         << "\n ---> The calorimeter is " << fNbOfLayers << " layers of:";
  for (G4int i = 1; i <= fNbOfAbsor; ++i) {
    G4cout << "\n \t" << std::setw(12) << fAbsorMaterial[i]->GetName() << ": " << std::setw(6)
           << G4BestUnit(fAbsorThickness[i], "Length");
  }
  G4cout << "\n-------------------------------------------------------------\n";

  G4cout << "\n" << fWorldMaterial << G4endl;
  for (G4int j = 1; j <= fNbOfAbsor; ++j) {
    G4cout << "\n" << fAbsorMaterial[j] << G4endl;
  }
  G4cout << "\n-------------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& material)
{
  // search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fWorldMaterial = pttoMaterial;
    if (fLogicWorld) {
      fLogicWorld->SetMaterial(fWorldMaterial);
      fLogicLayerFront->SetMaterial(fWorldMaterial);
      fLogicLayerBack->SetMaterial(fWorldMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfLayers(G4int ival)
{
  // set the number of Layers
  //
  if (ival < 2) {
    G4cout << "\n --->warning from SetfNbOfLayers: " << ival
           << " must be at least 2. Command refused" << G4endl;
    return;
  }
  fNbOfLayers = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNbOfAbsor(G4int ival)
{
  // set the number of Absorbers
  //
  if (ival < 1 || ival > (kMaxAbsor - 1)) {
    G4cout << "\n ---> warning from SetfNbOfAbsor: " << ival << " must be at least 1 and and most "
           << kMaxAbsor - 1 << ". Command refused" << G4endl;
    return;
  }
  fNbOfAbsor = ival;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorMaterial(G4int ival, const G4String& material)
{
  // search the material by its name
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorMaterial: absor number " << ival
           << " out of range. Command refused" << G4endl;
    return;
  }

  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(material);
  if (pttoMaterial) {
    fAbsorMaterial[ival] = pttoMaterial;
    if (fLogicAbsorFront[ival]) {
      fLogicAbsorFront[ival]->SetMaterial(pttoMaterial);
      fLogicAbsorBack[ival]->SetMaterial(pttoMaterial);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorThickness(G4int ival, G4double val)
{
  // change Absorber thickness
  //
  if (ival > fNbOfAbsor || ival <= 0) {
    G4cout << "\n --->warning from SetAbsorThickness: absor number " << ival
           << " out of range. Command refused" << G4endl;
    return;
  }
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetAbsorThickness: thickness " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fAbsorThickness[ival] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCalorSizeYZ(G4double val)
{
  // change the transverse size
  //
  if (val <= DBL_MIN) {
    G4cout << "\n --->warning from SetfCalorSizeYZ: thickness " << val
           << " out of range. Command refused" << G4endl;
    return;
  }
  fCalorSizeYZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AutoDelete.hh"
#include "G4GlobalMagFieldMessenger.hh"

void DetectorConstruction::ConstructSDandField()
{
  if (fFieldMessenger.Get() == nullptr) {
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    G4GlobalMagFieldMessenger* msg = new G4GlobalMagFieldMessenger(fieldValue);
    // msg->SetVerboseLevel(1);
    G4AutoDelete::Register(msg);
    fFieldMessenger.Put(msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
