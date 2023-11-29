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
/// \file scavenger/src/DetectorConstruction.cc
/// \brief Implementation of the scavenger::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "ScoreSpecies.hh"
#include "PrimaryKiller.hh"

namespace scavenger
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VPhysicalVolume *DetectorConstruction::Construct() {
  // Water is defined from NIST material database
  auto *man = G4NistManager::Instance();
  auto *water = man->FindOrBuildMaterial("G4_WATER");

  // World
  G4double world_sizeXYZ = 1. * km;
  auto solidWorld =
    new G4Box("World",
              0.5 * world_sizeXYZ,
              0.5 * world_sizeXYZ,
              0.5 * world_sizeXYZ);

  auto logicWorld =
    new G4LogicalVolume(solidWorld,
                        water,
                        "World");

  G4VPhysicalVolume *physWorld =
    new G4PVPlacement(nullptr,               //no rotation
                      G4ThreeVector(),       //its position at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      nullptr,               //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      true);                 //checking overlaps

  // Visualization attributes
  auto worldVisAtt = new G4VisAttributes(G4Colour(.5, 1.0, .5));
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void DetectorConstruction::ConstructSDandField() {
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  /**
   * declare World as a MultiFunctionalDetector scorer
   */
  auto mfDetector = new G4MultiFunctionalDetector("mfDetector");

  // Kill primary track after a chosen energy loss OR under a chosen
  // kinetic energy
  auto primaryKiller = new PrimaryKiller("PrimaryKiller");
  primaryKiller->SetMinLossEnergyLimit(500. * eV); // default value
  primaryKiller->SetMaxLossEnergyLimit(1. * keV); // default value
  mfDetector->RegisterPrimitive(primaryKiller);

  // Record Species scorer:
  //  - score number of species over time
  //  - score the total energy deposition
  //  - compute the radiochemical yields (G values)
  G4VPrimitiveScorer *gValues = new ScoreSpecies("Species");
  mfDetector->RegisterPrimitive(gValues);
  G4SDManager::GetSDMpointer()->AddNewDetector(mfDetector);
  SetSensitiveDetector("World", mfDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}