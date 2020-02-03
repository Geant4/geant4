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

#include "DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), fDiameter(6*nm), fLength(10*nm)
{
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct() {
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
    // Cleanup old geometry
    //
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    //
    G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");

    // Geometry parameters
    G4double worldSize = 150 * nanometer;

    G4VSolid* world = new G4Box("world",           // its name
            worldSize/2, worldSize/2, worldSize/2); // its size

    G4LogicalVolume* worldLog = new G4LogicalVolume(
            world,            // its solid
            water,            // its material
            "world");         // its name

    G4VisAttributes *vis = new G4VisAttributes(G4Colour(G4Colour::Blue()));
    worldLog->SetVisAttributes(vis);

    G4VPhysicalVolume* worldPhys = new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(),  // at (0,0,0)
            worldLog,         // its logical volume
            "world",          // its name
            0,                // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            true);            // checking overlaps forced to YES

    G4VSolid* target = new G4Tubs("target",
            0, fDiameter/2, fLength/2, 0, 360*deg); //its size

    G4LogicalVolume* targetLog = new G4LogicalVolume(
            target,
            water,
            "target");

    new G4PVPlacement(0,     // no rotation
            G4ThreeVector(), // at (0,0,0)
            targetLog,       // its logical volume
            "target",        // its name
            worldLog,        // its mother volume 
            false,           // no boolean operation
            0,               // copy number
            true);           // checking overlaps force to YES

    // A region is needed as splitting is performed in regions
    G4Region* aRegion = G4RegionStore::GetInstance()->FindOrCreateRegion("Target");
    targetLog->SetRegion(aRegion);
    aRegion->AddRootLogicalVolume(targetLog);

    G4VisAttributes *visTarget = new G4VisAttributes(G4Colour(G4Colour::Blue()));
    visTarget->SetForceSolid(true);
    visTarget->SetVisibility(true);
    targetLog->SetVisAttributes(visTarget);

    return worldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDiameter(G4double diameter) {
    fDiameter = diameter;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLength(G4double length) {
    fLength = length;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
    G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
