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
///  \file DetectorConstruction.cc
///  \brief Example of empty world volume.
//
//  Author: G.Hugo, 06 January 2023
//
// ***************************************************************************
//
///      DetectorConstruction
//
//  Example of empty world volume.
//
// ***************************************************************************

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SystemOfUnits.hh"


// ***************************************************************************
// Returns a 1-cm sided box filled with G4_Galactic.
// ***************************************************************************
G4VPhysicalVolume* DetectorConstruction::Construct() {

  
  G4Box* const worldSolid = new G4Box("World", 1.*CLHEP::cm, 1.*CLHEP::cm, 1.*CLHEP::cm);

  G4Material* const worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume* const worldLogicalVol = new G4LogicalVolume(worldSolid, worldMaterial, "World");
  // NB: G4LogicalVolumeStore owns all logical volumes.

  G4VPhysicalVolume* const worldPhysicalVol  = new G4PVPlacement(nullptr, 
                                                              G4ThreeVector(),
                                                              worldLogicalVol, 
                                                              "World", 
                                                              nullptr, 
                                                              false, 
                                                              0);
  // NB: G4PhysicalVolumeStore owns all physical volumes.

  return worldPhysicalVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
