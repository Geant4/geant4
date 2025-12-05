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

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Water is defined from NIST material database
  auto man = G4NistManager::Instance();
  auto water = man->FindOrBuildMaterial("G4_WATER");

  //
  // World
  //
  const G4double worldXYZ = 1 * m;

  auto solidWorld = new G4Box("World", 0.5 * worldXYZ, 0.5 * worldXYZ, 0.5 * worldXYZ);
  auto lvWorld = new G4LogicalVolume(solidWorld, water, "World");
  auto pvWorld =
    new G4PVPlacement(nullptr, G4ThreeVector(), lvWorld, "World", nullptr, false, 0, true);

  //
  // Cell
  ConstructCell(pvWorld);

  // return the world
  return pvWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructCell(G4VPhysicalVolume* pvWorld)
{
  const G4double cellRadius = 7 * um;
  const G4double nucleusRadius = 4 * um;

  const G4int nMitochondria = 100;
  const G4double mitoA = 0.55 * micrometer;
  const G4double mitoB = 0.25 * micrometer;
  const G4double mitoC = 0.90 * micrometer;

  auto solidCell = new G4Orb("Cell", cellRadius);
  auto lvCell = new G4LogicalVolume(
    solidCell, G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER"), "Cell");
  auto pvCell =
    new G4PVPlacement(nullptr, G4ThreeVector(), "Cell", lvCell, pvWorld, false, 0, true);

  auto solidNucleus = new G4Orb("Nucleus", nucleusRadius);
  auto lvNucleus = new G4LogicalVolume(
    solidNucleus, G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER"), "Nucleus");

  new G4PVPlacement(nullptr, G4ThreeVector(), "Nucleus", lvNucleus, pvCell, false, 0, true);

  auto solidMito = new G4Ellipsoid("Mitochondria", mitoA, mitoB, mitoC);
  auto lvMito = new G4LogicalVolume(
    solidMito, G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER"), "Mitochondria");

  for (auto i = 0; i < nMitochondria; ++i) {
    G4bool overlap = true;
    do {
      auto u = twopi * G4UniformRand();
      auto v = std::acos(2 * G4UniformRand() - 1);
      auto dr = G4UniformRand() * cellRadius;
      auto x = dr * std::cos(u) * std::sin(v);
      auto y = dr * std::sin(u) * std::sin(v);
      auto z = dr * std::cos(v);
      auto pos = G4ThreeVector(x, y, z);

      auto phi = G4UniformRand() * 2 * pi;
      auto psi = G4UniformRand() * 2 * pi;
      auto rot = new G4RotationMatrix();
      rot->rotateX(psi);
      rot->rotateY(phi);

      auto pvMito = new G4PVPlacement(rot, pos, "Mitochondria", lvMito, pvCell, false, i, false);

      overlap = pvMito->CheckOverlaps(1000, 0, false);
      if (overlap) {
        G4PhysicalVolumeStore::DeRegister(pvMito);
      }
    } while (overlap);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......