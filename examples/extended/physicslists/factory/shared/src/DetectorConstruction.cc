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
/// \file DetectorConstruction
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "ScreenSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();

  // Build materials
  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* csi = nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
       // There is no need to test if materials were built/found
       // as G4NistManager would issue an error otherwise
       // Try the code with "XYZ".      

  // Print all materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;         
  
  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  //
  // World
  //

  // The world dimensions
  G4double worldHxyz = 2.*m;
  
  // world volume
  G4Box* worldS = new G4Box("World", worldHxyz, worldHxyz, worldHxyz); 

  G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, air, "World");
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
            0, G4ThreeVector(), worldLV, "World", 0, false, 0, checkOverlaps);

  //
  // Box
  //

  // The box dimensions
  G4double boxHxy = 1.*m;
  G4double boxHz = 10.*cm;
  
  // box volume
  G4Box* boxS = new G4Box("World", boxHxy, boxHxy, boxHz); 

  G4LogicalVolume* boxLV = new G4LogicalVolume(boxS, csi, "Box");
                                   
  // The box position
  G4double posz = 0.*m;

  new G4PVPlacement(
        0, G4ThreeVector(0, 0, posz),
        boxLV, "Box", worldLV, false, 0, checkOverlaps);

  //
  // Scoring screen
  //
  // The screen dimensions
  G4double screenHxy = 1.999*m;
  G4double screenHz = 1.*mm;
  
  // Screen volume
  G4Box* screenS = new G4Box("World", screenHxy, screenHxy, screenHz); 
      
  G4LogicalVolume* screenLV = new G4LogicalVolume(screenS, air, "Screen");
                                   
  // The screen position
  posz += boxHz + screenHz;

  new G4PVPlacement(
        0, G4ThreeVector(0, 0, posz),
        screenLV, "Screen", worldLV, false, 0, checkOverlaps);

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  simpleBoxVisAtt->SetVisibility(true);
  boxLV->SetVisAttributes(simpleBoxVisAtt);
  screenLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //
  auto screenSD = new ScreenSD("ScreenSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(screenSD);
  SetSensitiveDetector("Screen", screenSD);

  // 
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
