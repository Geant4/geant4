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
/// \file GB01/src/GB01DetectorConstruction.cc
/// \brief Implementation of the GB01DetectorConstruction class
//
#include "GB01DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4NistManager.hh"

#include "GB01BOptrMultiParticleChangeCrossSection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB01DetectorConstruction::GB01DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB01DetectorConstruction::~GB01DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* GB01DetectorConstruction::Construct()
{
  G4Material*   worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4Material* defaultMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_lN2");

  G4VSolid* solidWorld = new G4Box("World", 10*m, 10*m, 10*m );
  
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,       // its solid
                                                    worldMaterial,    // its material
                                                    "World");         // its name
  
  G4PVPlacement* physiWorld = new G4PVPlacement(0,                    // no rotation
                                                G4ThreeVector(),      // at (0,0,0)
                                                logicWorld,           // its logical volume
                                                "World",              // its name
                                                0,                    // its mother volume
                                                false,                // no boolean operation
                                                0);                   // copy number
  
  // -----------------------------------
  // -- volume where biasing is applied:
  // -----------------------------------
  G4double halfZ = 10*cm;
  G4VSolid* solidTest = new G4Box("test.solid", 1*m, 1*m, halfZ );
  
  G4LogicalVolume* logicTest = new G4LogicalVolume(solidTest,        // its solid
                                                   defaultMaterial,  // its material
                                                   "test.logical");  // its name

  new G4PVPlacement(0,                                               // no rotation
                    G4ThreeVector(0,0, halfZ),                       // volume entrance at (0,0,0)
                    logicTest,                                       // its logical volume
                    "test.phys",                                     // its name
                    logicWorld,                                      // its mother  volume
                    false,                                           // no boolean operation
                    0);                                              // copy number
  
  
  return physiWorld;
}

void GB01DetectorConstruction::ConstructSDandField()
{
  // -- Fetch volume for biasing:
  G4LogicalVolume* logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume("test.logical");
  
  // ----------------------------------------------
  // -- operator creation and attachment to volume:
  // ----------------------------------------------
  GB01BOptrMultiParticleChangeCrossSection* testMany = 
    new GB01BOptrMultiParticleChangeCrossSection();
  testMany->AddParticle("gamma");
  testMany->AddParticle("neutron");
  testMany->AttachTo(logicTest);
  G4cout << " Attaching biasing operator " << testMany->GetName()
         << " to logical volume " << logicTest->GetName()
         << G4endl;
}
