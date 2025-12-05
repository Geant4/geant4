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
/// \brief Implementation of the RadioBio::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"

#include "DetectorMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

#include <sstream>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  // Set the default box material
  SetMaterial("G4_WATER");

  // Create the messenger
  fDetectorMessenger = new DetectorMessenger(this);

  // Voxelize the detecctor with default size
  VoxelizedSensitiveDetector::CreateInstance(this, 0.1 * m, 1 * m, 1 * m);
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

void DetectorConstruction::ConstructSDandField()
{
  VoxelizedSensitiveDetector::GetInstance()->ConstructSD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Materials
  G4bool isotopes = false;
  G4Material* airNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);

  // Simulation world
  const G4double worldX = 400.0 * cm;
  const G4double worldY = 400.0 * cm;
  const G4double worldZ = 400.0 * cm;

  G4Box* sWorld = new G4Box("TreatmentRoom", worldX, worldY, worldZ);

  G4LogicalVolume* lWorld = new G4LogicalVolume(sWorld, airNist, "logicWorld", 0, 0, 0);

  pWorld = new G4PVPlacement(0, G4ThreeVector(), "physicsWorld", lWorld, 0, false, 0);

  // Create a visual attribute for the world
  G4VisAttributes* worldVisAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));  // Blue color
  worldVisAttributes->SetVisibility(true);  // Ensure visibility
  worldVisAttributes->SetForceWireframe(true);  // Optional: make it wireframe for clarity
  lWorld->SetVisAttributes(worldVisAttributes);

  G4Box* sBox = new G4Box("Container", fBoxSizeX / 2, fBoxSizeY / 2, fBoxSizeZ / 2);

  fLBox = new G4LogicalVolume(sBox, fMaterial, fMaterial->GetName());

  fPBox = new G4PVPlacement(0, G4ThreeVector(), fLBox, fMaterial->GetName(), lWorld, false, 0);

  // Initialize pointer to world for voxelization
  VoxelizedSensitiveDetector::GetInstance()->InitializeWorldPtr(fPBox);

  // Create Voxelized Geometry
  VoxelizedSensitiveDetector::GetInstance()->Construct();

  // Always return the root volume
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box dimensions are: " << G4endl
         << "x: " << G4BestUnit(fBoxSizeX, "Length") << G4endl
         << "y: " << G4BestUnit(fBoxSizeY, "Length") << G4endl
         << "z: " << G4BestUnit(fBoxSizeZ, "Length") << G4endl;

  G4cout << "And its volume therefore is: "
         << fPBox->GetLogicalVolume()->GetSolid()->GetCubicVolume() / m3 << " m3" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    if (fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if (fLBox) {
        fLBox->SetMaterial(pttoMaterial);
      }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    std::stringstream sstr;
    sstr << "material " << +materialChoice << " does not exist, keeping material "
         << fMaterial->GetName();

    G4Exception("DetectorConstruction::SetMaterial", "NoWorldMat", JustWarning, sstr.str().c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  SetSizeX(value);
  SetSizeY(value);
  SetSizeZ(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4ThreeVector size)
{
  SetSizeX(size.getX());
  SetSizeY(size.getY());
  SetSizeZ(size.getZ());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double sizeX)
{
  fBoxSizeX = sizeX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeY(G4double sizeY)
{
  fBoxSizeY = sizeY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeZ(G4double sizeZ)
{
  fBoxSizeZ = sizeZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio