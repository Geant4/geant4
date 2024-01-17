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
//
/// \file radiobiology/src/DetectorConstruction.cc
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

  G4Box* sBox = new G4Box("Container",  // its name
                          fBoxSizeX / 2, fBoxSizeY / 2, fBoxSizeZ / 2);  // its dimensions

  fLBox = new G4LogicalVolume(sBox,  // its shape
                              fMaterial,  // its material
                              fMaterial->GetName());  // its name

  fPBox = new G4PVPlacement(0,  // no rotation
                            G4ThreeVector(),  // at (0,0,0)
                            fLBox,  // its logical volume
                            fMaterial->GetName(),  // its name
                            0,  // its mother  volume
                            false,  // no boolean operation
                            0);  // copy number

  // Parameters for the world volume can be printed
  // PrintParameters();

  // Initialize pointer to world for voxelization
  VoxelizedSensitiveDetector::GetInstance()->InitializeWorldPtr(fPBox);

  // Create Voxelized Geometry
  VoxelizedSensitiveDetector::GetInstance()->Construct();

  // Always return the root volume
  return fPBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box dimensions are: " << G4endl << "x: " << G4BestUnit(fBoxSizeX, "Length")
         << G4endl << "y: " << G4BestUnit(fBoxSizeY, "Length") << G4endl
         << "z: " << G4BestUnit(fBoxSizeZ, "Length") << G4endl;

  G4cout << "And its volume therefore is: "
         << fPBox->GetLogicalVolume()->GetSolid()->GetCubicVolume() / m3 << " m3" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // Search the material by its name
  G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    if (fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if (fLBox) {
        fLBox->SetMaterial(pttoMaterial);
      }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
  else {
    // Warning the user this material does not exist
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
