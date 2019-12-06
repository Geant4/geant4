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
/// \file GB07/src/GB07DetectorConstruction.cc
/// \brief Implementation of the GB07DetectorConstruction class
//
#include "GB07DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"

#include "G4NistManager.hh"

#include "GB07BOptrLeadingParticle.hh"

#include "GB07SD.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07DetectorConstruction::GB07DetectorConstruction()
  : G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB07DetectorConstruction::~GB07DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GB07DetectorConstruction::Construct()
{
  G4Material*   worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4Material* defaultMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
  
  G4VSolid* solidWorld = new G4Box("World", 10*m, 10*m, 10*m );
  
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,         // its solid
                                                    worldMaterial,      // its material
                                                    "World");           // its name
  
  G4PVPlacement* physiWorld = new G4PVPlacement(0,                      // no rotation
                                                G4ThreeVector(),        // at (0,0,0)
                                                logicWorld,             // its logical volume
                                                "World",                // its name
                                                0,                      // its mother  volume
                                                false,                  // no boolean operation
                                                0);                     // copy number
  
  // -----------------------------------
  // -- volume where biasing is applied:
  // -----------------------------------
  G4double halfZ = 1*m;
  G4VSolid* solidTest = new G4Box("test.solid", 1*m, 1*m, halfZ );
  
  G4LogicalVolume* logicTest = new G4LogicalVolume(solidTest,           // its solid
                                                   defaultMaterial,     // its material
                                                   "test.logical");     // its name
  
  new G4PVPlacement(0,                                                  // no rotation
                    G4ThreeVector(0,0, halfZ),                          // put entrance at (0,0,0)
                    logicTest,                                          // its logical volume
                    "test.physical",                                    // its name
                    logicWorld,                                         // its mother  volume
                    false,                                              // no boolean operation
                    0);                                                 // copy number
  
  
  // ------------------------------------------------------------
  // -- volume to tally exiting particles, put next to above one:
  // ------------------------------------------------------------
  G4double halfZtally = 0.5*mm;
  G4VSolid* solidTally = new G4Box("tally.solid", 1*m, 1*m, halfZtally );
  G4LogicalVolume* logicTally = new G4LogicalVolume(solidTally,         // its solid
                                                    worldMaterial,      // its material
                                                    "tally.logical");   // its name
   
  new G4PVPlacement(0,                                                  // no rotation
                    G4ThreeVector(0,0, 2*halfZ + halfZtally),           // put next to test.phys
                    logicTally,                                         // its logical volume
                    "tally.physical",                                   // its name
                    logicWorld,                                         // its mother volume
                    false,                                              // no boolean operation
                    0);                                                 // copy number
 
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB07DetectorConstruction::ConstructSDandField()
{
  // -- Fetch volume for biasing:
  G4LogicalVolume* logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("test.logical");
  
  // ----------------------------------------------
  // -- operator creation and attachment to volume:
  // ----------------------------------------------
  GB07BOptrLeadingParticle* testMany =  new GB07BOptrLeadingParticle();
  testMany->AttachTo(logicTest);
  G4cout << " Attaching biasing operator " << testMany->GetName()
         << " to logical volume " << logicTest->GetName()
         << G4endl;
  
  // ---------------------------------------------------------------------------------
  // -- Attach sensitive detector to print information on particle exiting the shield:
  // ---------------------------------------------------------------------------------
  // -- create and register sensitive detector code module:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* sd = new GB07SD("Tally");
  SDman->AddNewDetector(sd);
  // -- Fetch volume for sensitivity and attach sensitive module to it:
  G4LogicalVolume* logicSD =
    G4LogicalVolumeStore::GetInstance()->GetVolume("tally.logical");
  logicSD->SetSensitiveDetector(sd);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
