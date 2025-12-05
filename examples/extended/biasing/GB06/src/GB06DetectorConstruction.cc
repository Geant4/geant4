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
/// \file GB06DetectorConstruction.cc
/// \brief Implementation of the GB06DetectorConstruction class

#include "GB06DetectorConstruction.hh"

#include "GB06SD.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06DetectorConstruction::GB06DetectorConstruction() : G4VUserDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06DetectorConstruction::~GB06DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GB06DetectorConstruction::Construct()
{
  G4Material* worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4Material* concreteMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");

  G4VSolid* solidWorld = new G4Box("World.solid", 10 * m, 10 * m, 10 * m);

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        worldMaterial,  // its material
                                        "World.logical");  // its name

  auto physiWorld = new G4PVPlacement(nullptr,  // no rotation
                                      G4ThreeVector(),  // at (0,0,0)
                                      logicWorld,  // its logical volume
                                      "World.physical",  // its name
                                      nullptr,  // its mother volume
                                      false,  // no bool. operation
                                      0);  // copy number

  // ----------------------------------------------------
  // -- volume of shield, made of concrete, in one block:
  // ----------------------------------------------------
  G4double halfXY = 1.5 * m;
  G4double halfZ = 2.5 * m;
  G4VSolid* solidShield = new G4Box("shield.solid", halfXY, halfXY, halfZ);

  auto logicTest = new G4LogicalVolume(solidShield,  // its solid
                                       concreteMaterial,  // its material
                                       "shield.logical");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(0, 0, halfZ),  // volume entrance is set at (0,0,0)
                    logicTest,  // its logical volume
                    "shield.physical",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0);  // copy number

  // ------------------------------------------------------------
  // -- dummy volume to display exiting neutron flux information:
  // ------------------------------------------------------------
  G4double halfz = 1 * cm;
  G4VSolid* solidMeasurement = new G4Box("meas.solid", halfXY, halfXY, halfz);

  logicMeasurement = new G4LogicalVolume(solidMeasurement,  // its solid
                                         worldMaterial,  // its material
                                         "meas.logical");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(0, 0, 2 * halfZ + halfz),  // entrance set after shield
                    logicMeasurement,  // its logical volume
                    "meas.physical",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0);  // copy number

  // -- world volume pointer returned:
  return physiWorld;
}

void GB06DetectorConstruction::ConstructSDandField()
{
  // ---------------------------------------------------------------------------------
  // -- Attach sensitive detector to print information on particle exiting the shield:
  // ---------------------------------------------------------------------------------
  // -- create and register sensitive detector code module:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* sd = new GB06SD("Measurer");
  SDman->AddNewDetector(sd);
  // -- Fetch volume for sensitivity and attach sensitive module to it:
  logicMeasurement->SetSensitiveDetector(sd);
}
