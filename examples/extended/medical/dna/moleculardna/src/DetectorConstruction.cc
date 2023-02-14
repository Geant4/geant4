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
/// \file DetectorConstruction.cc
/// \brief Detector Construction for a DNA molecule in a cell

#include "DetectorConstruction.hh"
#include "DNAGeometry.hh"
#include "DetectorMessenger.hh"
#include "OctreeNode.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
  , fpDNAGeometry(new DNAGeometry())
  , fpDetectorMessenger(new DetectorMessenger(this))
{
  G4bool useParallelPhysicsWorld = false;
  if (useParallelPhysicsWorld) {
    RegisterParallelWorld(fpDNAGeometry->GetDNAWorld());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() { delete fpDetectorMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  BuildMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::BuildMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  fpWorld            = man->FindOrBuildMaterial("G4_Galactic");
  fpWater            = man->FindOrBuildMaterial("G4_WATER");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  //////////////////////////////////////////////////////////////////////////////
  // Brief :: Build the detector. This real world for only chemistry stage only.
  // An other parallel world is built in DNAWorld for only physics We either
  // build a DNA region inside a cell if fCellRadius is se in macro file, or we
  // just run the experiment in the World if it isn't. If we run in a cell, the
  // world material is a vacuum, and the cell is water. Otherwise everything is
  // water.
  //////////////////////////////////////////////////////////////////////////////

  G4LogicalVolume* pDnaRegionLogical;

  G4double dnasx = fCellRadius.getX();
  G4double dnasy = fCellRadius.getY();
  G4double dnasz = fCellRadius.getZ();
  if((dnasx <= 0) || (dnasy <= 0) || (dnasz <= 0))
  {
    // Experiment in the world only
    G4cout
      << "No cell radius has been specified (or the radii are less than 0)."
      << G4endl;
    G4cout << "This will fill the entire world with water" << G4endl;

    G4double wsize     = 1.01 * fWorldSideLength / 2.;
    auto worldPhysical = new G4Box("WorldPhysical", wsize, wsize, wsize);
    auto worldLogical  = new G4LogicalVolume(worldPhysical, fpWater, "world");
    fWorld = new G4PVPlacement(nullptr, G4ThreeVector(0), worldLogical, "world",
                               nullptr, false, 0);
    pDnaRegionLogical = worldLogical;
  }
  else
  {
    // experiment inside a cell.
    G4double wsize     = 1.01 * fWorldSideLength / 2.;
    auto worldPhysical = new G4Box("WorldPhysical", wsize, wsize, wsize);
    auto worldLogical  = new G4LogicalVolume(worldPhysical, fpWorld, "world");

    fWorld = new G4PVPlacement(nullptr, G4ThreeVector(0), worldLogical, "world",
                               nullptr, false, 0);
    auto dnaPhysical = new G4Ellipsoid("CellPhysical", dnasx, dnasy, dnasz);
    auto dnaLogical  = new G4LogicalVolume(dnaPhysical, fpWater, "CellLogical");

    new G4PVPlacement(nullptr, G4ThreeVector(0), "CellVolume", dnaLogical,
                      fWorld, false, 0, false);
    pDnaRegionLogical = dnaLogical;
  }
  // we build DNA geometry for physical stage
  this->GetDNAGeometry()->BuildDNA(pDnaRegionLogical);
  return fWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
