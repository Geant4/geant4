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
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():
fpDefaultMaterial(0),fpWaterMaterial(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();

  nistManager->FindOrBuildMaterial("G4_WATER", false);
  fpWaterMaterial = G4Material::GetMaterial("G4_WATER");

  nistManager->FindOrBuildMaterial("G4_AIR", false);
  fpDefaultMaterial = G4Material::GetMaterial("G4_AIR");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // World
  G4double worldSize  = 2*micrometer;

  G4Box* solidWorld = new G4Box("World",
                                worldSize/2.,
                                worldSize/2.,
                                worldSize/2.);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                    fpDefaultMaterial,
                                                    "World");
  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,
                                                    G4ThreeVector(),
                                                    "World",
                                                    logicWorld,
                                                    0,
                                                    false,
                                                    0);
  // Target Box
  // from Francis et al. Comput. Meth. Programs. Biomed. 2011 101(3)
  G4double targetSizeXY=1*micrometer;
  G4double targetSizeZ=0.5*micrometer;

  G4Box* solidTarget
  = new G4Box("Target", targetSizeXY/2.,targetSizeXY/2.,targetSizeZ/2.);
  G4LogicalVolume* logicTarget
  = new G4LogicalVolume(solidTarget, fpWaterMaterial, "Target");
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    "Target",
                    logicTarget,
                    physiWorld,
                    false,
                    0);

  // Visualization attributes
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);

  G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  targetVisAtt->SetVisibility(true);
  targetVisAtt->SetForceAuxEdgeVisible(true);
  logicTarget->SetVisAttributes(targetVisAtt);

  //always return the world volume
  return physiWorld;
}

