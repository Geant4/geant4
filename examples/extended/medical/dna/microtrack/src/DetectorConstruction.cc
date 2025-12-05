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

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4GeometryManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  // Create commands for interactive definition of the detector
  fDetectorMessenger = std::make_unique<DetectorMessenger>(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define the materials
  DefineMaterials();

  // Define volumes
  G4VPhysicalVolume* World = DefineWorld();
  DefineSD(World);  // SD is placed inside the world volume
  return World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Taking water from the NIST database
  G4NistManager* nist = G4NistManager::Instance();
  if (!fMat) {
    fMat = nist->FindOrBuildMaterial("G4_WATER");
  }

  // Print the material information
  G4cout << "Material: " << fMat->GetName() << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckConsistency()
{
  if (!(fMaxRange > 0.)) {
    G4ExceptionDescription msg;
    msg << "fMaxRange must be > 0.\n"
        << "Please consider using the formula by Tabata for the calculation "
        << "of the maximum range of secondary electrons:\n"
        << "\t https://doi.org/10.1016/0029-554X(72)90463-6" << G4endl
        << "Note: This parameter affects the world volume size only along "
        << "the z-axis.\n"
        << "Sizes along x- and y-axes are set by fHitSelRegXY and the "
        << "site radius." << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0001",
                FatalException, msg);
  }
  if (!(fHitSelRegZ > 0.)) {
    G4ExceptionDescription msg;
    msg << "fHitSelRegZ must be > 0.\n"
        << "Otherwise, no hits are eligible to set randomly the site "
        << "position using the weighted method.\n"
        << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0002",
                FatalException, msg);
  }
  if (!(fHitSelRegXY > 0.)) {
    G4ExceptionDescription msg;
    msg << "fHitSelRegXY must be > 0.\n"
        << "Otherwise, no hits are eligible to set randomly the site "
        << "position using the weighted method." << G4endl
        << "Note: This parameter, together with the site radius, sets the "
        << "world volume size along the xy-axes." << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0003",
                FatalException, msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineWorld()
{
  // Ensure geometry consistency
  CheckConsistency();

  // Define the world volume

  G4double WorldZ = fHitSelRegZ + 4. * fSiteRadius + 2. * fMaxRange;
  G4double WorldXY = fHitSelRegXY + 4. * fSiteRadius;

  auto* sWorld = new G4Box("World", WorldXY / 2., WorldXY / 2., WorldZ / 2.);
  auto* lWorld = new G4LogicalVolume(sWorld, fMat, "World", 0, 0, 0);
  G4VPhysicalVolume* World =
    new G4PVPlacement(0, G4ThreeVector(), lWorld, "World", nullptr, false, 0);

  // Print the world volume information
  PrintParameters(World);

  return World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
DetectorConstruction::DefineSD(G4VPhysicalVolume* mother) const
{
  // Define the sensitive detector volume

  G4double SensDetZ = fHitSelRegZ + 4. * fSiteRadius;
  G4double SensDetXY = fHitSelRegXY + 4. * fSiteRadius;

  auto* sSDbox =
    new G4Box("SDbox", SensDetXY / 2., SensDetXY / 2., SensDetZ / 2.);
  auto* lSDbox = new G4LogicalVolume(sSDbox, fMat, "SDbox", 0, 0, 0);
  G4VPhysicalVolume* SDbox = new G4PVPlacement(
    0, G4ThreeVector(), lSDbox, "SDbox", mother->GetLogicalVolume(), false, 0);
  // Print the sensitive detector volume information
  PrintParameters(SDbox);

  return SDbox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive Detector
  G4String SDname = "SDbox";

  auto* aSD = new TrackerSD(SDname, "TrackerHitColl");

  // Register the SD with the SD manager
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  SDMan->AddNewDetector(aSD);

  // Set the SD to the logical volume
  SetSensitiveDetector("SDbox", aSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& name)
{
  auto* nist = G4NistManager::Instance();

  G4Material* mat = nist->FindOrBuildMaterial(name);

  if (!mat) {
    G4ExceptionDescription ed;
    ed << "Material '" << name << "' not found. Falling back to G4_WATER.";
    G4Exception("DetectorConstruction::SetMaterial", "MTK001", JustWarning, ed);
    mat = nist->FindOrBuildMaterial("G4_WATER");
  }

  if (fMat != mat) {
    fMat = mat;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters(G4VPhysicalVolume* physVol) const
{
  if (!physVol) {
    G4cerr << "Error: Physical volume is null!" << G4endl;
    return;
  }

  G4ThreeVector pos = physVol->GetObjectTranslation();
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  G4VSolid* solid = logVol->GetSolid();
  G4Material* material = logVol->GetMaterial();

  G4cout << "\n================ Volume Parameters ================" << G4endl;
  G4cout << "Physical Volume Name: " << physVol->GetName() << G4endl;
  G4cout << "Position: " << pos / mm << " mm" << G4endl;
  G4cout << "Material: " << material->GetName() << G4endl;

  // Detect shape and print dimensions accordingly
  if (auto box = dynamic_cast<G4Box*>(solid)) {
    G4cout << "Shape: Box" << G4endl;
    G4cout << "Dimensions (full): " << 2 * box->GetXHalfLength() / mm << " x "
           << 2 * box->GetYHalfLength() / mm << " x "
           << 2 * box->GetZHalfLength() / mm << " mm" << G4endl;
  }
  else {
    G4cout << "Shape: Unknown or not handled yet." << G4endl;
  }

  G4cout << "===================================================\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
