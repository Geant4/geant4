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
/// \file ExG4DetectorConstruction01.cc
/// \brief Implementation of the ExG4DetectorConstruction01 class

#include "ExG4DetectorConstruction01.hh"
#include "ExG4DetectorConstruction01Messenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction01::ExG4DetectorConstruction01(
                              const G4String& materialName,
                              G4double hx, G4double hy, G4double hz)
 : G4VUserDetectorConstruction(),
   fMessenger(this),
   fMaterialName(materialName),
   fDimensions(hx, hy, hz),
   fWorldVolume(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4DetectorConstruction01::~ExG4DetectorConstruction01()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExG4DetectorConstruction01::Construct()
{
  // Define materials via NIST manager
  //
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* material
    = nistManager->FindOrBuildMaterial(fMaterialName);
 
  // World
  //
  G4Box* sWorld 
    = new G4Box("World",                        //name
                 fDimensions.x(),               //dimensions (half-lentghs)
                 fDimensions.y(), 
                 fDimensions.z());

  fWorldVolume 
    = new G4LogicalVolume(sWorld,               //shape
                          material,             //material
                          "World");             //name

  G4VPhysicalVolume* pWorld 
    = new G4PVPlacement(0,                      //no rotation
                        G4ThreeVector(),        //at (0,0,0)
                        fWorldVolume,           //logical volume
                        "World",                //name
                        0,                      //mother  volume
                        false,                  //no boolean operation
                        0);                     //copy number

  //always return the root volume
  //
  return pWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExG4DetectorConstruction01::SetMaterial(const G4String& materialName)
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
 
void ExG4DetectorConstruction01::SetDimensions(
                                   G4double hx, G4double hy, G4double hz) 
{
/// Set world dimension (in half lengths).
/// This setting has effect only if called in PreInit> phase

  fDimensions = G4ThreeVector(hx, hy, hz);
}  
                                     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
