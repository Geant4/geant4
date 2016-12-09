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
// $Id$
//
/// \file ExG4DetectorConstruction02.cc
/// \brief Implementation of the ExG4DetectorConstruction02 class

#include "ExG4DetectorConstruction02.hh"
#include "ExG4DetectorConstruction02Messenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction02::ExG4DetectorConstruction02(
                              const G4String& boxMaterialName,
                              G4double boxHx, G4double boxHy, G4double boxHz,
                              const G4String& worldMaterialName,
                              G4double worldSizeFactor)
 : G4VUserDetectorConstruction(),
   fMessenger(this),
   fBoxMaterialName(boxMaterialName),
   fWorldMaterialName(worldMaterialName),
   fBoxDimensions(boxHx*2, boxHy*2, boxHz*2),
   fWorldSizeFactor(worldSizeFactor),
   fBoxVolume(0),
   fWorldVolume(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction02::~ExG4DetectorConstruction02()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExG4DetectorConstruction02::Construct()
{
  // Define materials via NIST manager
  //
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* worldMaterial
    = nistManager->FindOrBuildMaterial(fWorldMaterialName);

  G4Material* boxMaterial
    = nistManager->FindOrBuildMaterial(fBoxMaterialName);
 
  // Geometry parameters
  //
  G4ThreeVector worldDimensions = fBoxDimensions * fWorldSizeFactor;
  
  
  // World
  //
  G4Box* sWorld 
    = new G4Box("World",                        //name
                 worldDimensions.x(),           //dimensions (half-lentghs)
                 worldDimensions.y(), 
                 worldDimensions.z());

  fWorldVolume 
    = new G4LogicalVolume(sWorld,               //shape
                          worldMaterial,        //material
                          "World");             //name

  G4VPhysicalVolume* pWorld 
    = new G4PVPlacement(0,                      //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                        fWorldVolume,           //logical volume
                        "World",                //name
                        0,                      //mother  volume
                        false,                  //no boolean operation
                        0);                     //copy number
                                                
  // Box
  //                           
  G4Box* sBox 
    = new G4Box("Box",                          //its name
                 fBoxDimensions.x(),            //dimensions (half-lengths)
                 fBoxDimensions.y(), 
                 fBoxDimensions.z());
                   
  fBoxVolume 
    = new G4LogicalVolume(sBox,                 //its shape
                          boxMaterial,          //its material
                          "Box");               //its name

  new G4PVPlacement(0,                          //no rotation
                    G4ThreeVector(),            //at (0,0,0)
                    fBoxVolume,                 //its logical volume                           
                    "Box",                      //its name
                    fWorldVolume,               //its mother  volume
                    false,                      //no boolean operation
                    0);                         //copy number

  //always return the root volume
  //
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExG4DetectorConstruction02::SetBoxMaterial(const G4String& materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* newMaterial
    = nistManager->FindOrBuildMaterial(materialName);

  if ( ! newMaterial ) {
    G4cerr << "Material " << materialName << " not found." << G4endl;
    G4cerr << "The box material was not changed." << G4endl;
    return;
  }  
   
  if ( fBoxVolume ) fBoxVolume->SetMaterial(newMaterial);
  G4cout << "Material of box changed to " << materialName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExG4DetectorConstruction02::SetWorldMaterial(const G4String& materialName)
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* newMaterial
    = nistManager->FindOrBuildMaterial(materialName);

  if ( ! newMaterial ) {
    G4cerr << "Material " << materialName << " not found." << G4endl;
    G4cerr << "The box material was not changed." << G4endl;
    return;
  }  
   
  if ( fWorldVolume ) fWorldVolume->SetMaterial(newMaterial);
  G4cout << "Material of box changed to " << materialName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExG4DetectorConstruction02::SetBoxDimensions(
                                   G4double hx, G4double hy, G4double hz) 
{
/// Set box dimension (in half lengths).
/// This setting has effect only if called in PreInit> phase

  fBoxDimensions = G4ThreeVector(hx, hy, hz);
}  
                                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExG4DetectorConstruction02::SetWorldSizeFactor(G4double factor) 
{
/// Set the multiplication factor from box dimensions to world dimensions.
/// This setting has effect only if called in PreInit> phase

  fWorldSizeFactor = factor;
}  
                                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
