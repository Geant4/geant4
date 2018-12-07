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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "G4ProductionCuts.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Orb.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material *
DetectorConstruction::OtherMaterial(G4String materialName)
{
  G4Material * material(0);

  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  material = man->FindOrBuildMaterial(materialName);

  // If one wishes to test other density value for water material,
  // one should use instead:
  // G4Material * H2O =
  //  man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED",
  //      "G4_WATER",
  //      1000*g/cm/cm/cm);
  // Note: any string for "G4_WATER_MODIFIED" parameter is accepted
  // and "G4_WATER" parameter should not be changed

  return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  G4Material *water = OtherMaterial("G4_WATER");

  // WORLD VOLUME = an G4ORB FULL OF LIQUID WATER

  double worldRadius = 0.3*micrometer;

  G4Orb* solidWorld = new G4Orb("World", worldRadius);

  G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,  //its solid
          water,  //its material
          "World");  //its name

  G4VPhysicalVolume* physiWorld =
      new G4PVPlacement(0,  //no rotation
          G4ThreeVector(),  //at (0,0,0)
          "World",  //its name
          logicWorld,  //its logical volume
          0,  //its mother  volume
          false,  //no boolean operation
          0);  //copy number

  // Visualization attributes = semitransparent navy blue
  G4VisAttributes* worldVisAtt =
      new G4VisAttributes(G4Colour(0.0,0.4,0.8,0.5));
  worldVisAtt->SetForceSolid(true);
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);

  return physiWorld;
}
