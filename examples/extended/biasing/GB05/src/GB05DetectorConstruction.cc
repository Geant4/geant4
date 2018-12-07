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
/// \file GB05DetectorConstruction.cc
/// \brief Implementation of the GB05DetectorConstruction class

#include "GB05DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4NistManager.hh"

#include "GB05BOptrSplitAndKillByCrossSection.hh"

#include "GB05SD.hh"
#include "G4SDManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB05DetectorConstruction::GB05DetectorConstruction()
 : G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB05DetectorConstruction::~GB05DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GB05DetectorConstruction::Construct()
{
  // -- Collect a few materials from NIST database:
  auto nistManager = G4NistManager::Instance();
  G4Material*   worldMaterial = nistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* defaultMaterial = nistManager->FindOrBuildMaterial("G4_CONCRETE");


  G4VSolid* solidWorld = new G4Box("World", 10*m, 10*m, 10*m );
  
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,    // its solid
                                                    worldMaterial, // its material
                                                    "World");      // its name
  
  G4PVPlacement* physiWorld = new G4PVPlacement(0,                 // no rotation
                                                G4ThreeVector(),   // at (0,0,0)
                                                logicWorld,        // its logical volume
                                                "World",           // its name
                                                0,                 // its mother volume
                                                false,             // no boolean operation
                                                0);                // copy number
  
  // -----------------------------------
  // -- volume where biasing is applied:
  // -----------------------------------
  G4double halfXY = 5.0*m;
  G4double halfZ  = 1.0*m;
  G4VSolid* solidShield = new G4Box("shield.solid", halfXY, halfXY, halfZ );
  
  G4LogicalVolume* logicShield = new G4LogicalVolume( solidShield,       // its solid
                                                      defaultMaterial,   // its material
                                                      "shield.logical"); // its name

  new G4PVPlacement(0,                                      // no rotation
                    G4ThreeVector(0,0, halfZ),              // volume entrance at (0,0,0)
                    logicShield,                            // its logical volume
                    "shield.phys",                          // its name
                    logicWorld,                             // its mother  volume
                    false,                                  // no boolean operation
                    0);                                     // copy number
  
  // --------------------------------------------------------------------
  // -- dummy volume to print information on particle exiting the shield:
  // --------------------------------------------------------------------
  G4double halfz = 1*cm;
  G4VSolid*        solidMeasurement = new G4Box("meas.solid", halfXY, halfXY, halfz );
  
  G4LogicalVolume* logicMeasurement = new G4LogicalVolume(solidMeasurement,// its solid
                                                          worldMaterial,   // its material
                                                          "meas.logical"); // its name
  
  new G4PVPlacement(0,                                      // no rotation
                    G4ThreeVector(0,0, 2*halfZ + halfz),    // volume entrance at (0,0,0)
                    logicMeasurement,                       // its logical volume
                    "meas.phys",                            // its name
                    logicWorld,                             // its mother  volume
                    false,                                  // no boolean operation
                    0);                                     // copy number
  
  
  return physiWorld;
}


void GB05DetectorConstruction::ConstructSDandField()
{
  // -- Fetch volume for biasing:
  G4LogicalVolume* logicShield =
    G4LogicalVolumeStore::GetInstance()->GetVolume("shield.logical");
  
  // -------------------------------------------------------------
  // -- operator creation, configuration and attachment to volume:
  // -------------------------------------------------------------
  GB05BOptrSplitAndKillByCrossSection* biasingOperator = 
    new GB05BOptrSplitAndKillByCrossSection("neutron");
  // -- Now, we declare to our biasing operator all the processes we
  // -- their disapperance effect on neutrons to be counterbalanced
  // -- by the splitting by cross-section :
  biasingOperator->AddProcessToEquipoise("Decay");
  biasingOperator->AddProcessToEquipoise("nCapture");
  biasingOperator->AddProcessToEquipoise("neutronInelastic");

  biasingOperator->AttachTo(logicShield);
  
  G4cout << " Attaching biasing operator " << biasingOperator->GetName()
         << " to logical volume " << biasingOperator->GetName()
         << G4endl;

  // ------------------------------------------------------------------------------------
  // -- Attach a sensitive detector to print information on particles exiting the shield:
  // ------------------------------------------------------------------------------------
  // -- create and register sensitive detector code module:
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4VSensitiveDetector* sd = new GB05SD("Scorer");
  SDman->AddNewDetector(sd);
  // -- Fetch volume for sensitivity and attach sensitive module to it:
  G4LogicalVolume* logicSD =
    G4LogicalVolumeStore::GetInstance()->GetVolume("meas.logical");
  logicSD->SetSensitiveDetector(sd);
  
}
