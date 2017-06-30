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
// $ID$
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fpWaterMaterial(0),fpRegion(0)
{}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()

{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 

  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");

  // Default materials in setup.
  fpWaterMaterial = H2O;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  // WORLD VOLUME
  
  G4double worldSizeX  = 1*mm;
  G4double worldSizeY  = 1*mm;
  G4double worldSizeZ  = 1*mm;

  G4VSolid* solidWorld = new G4Box("World",         //its name
                          worldSizeX/2,
                          worldSizeY/2,
                          worldSizeZ/2);  //its size

  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,  //its solid
                                    fpWaterMaterial,  //its material
                                    "World");    //its name

  G4VPhysicalVolume* physiWorld = new G4PVPlacement(0,      //no rotation
                                  G4ThreeVector(),  //at (0,0,0)
                                  "World",    //its name
                                  logicWorld,    //its logical volume
                                  0,      //its mother  volume
                                  false,      //no boolean operation
                                  0);      //copy number

  G4double TargetSizeZ =  worldSizeZ*0.05; 

  G4Box* targetSolid = new G4Box("Target",               //its name
                                 worldSizeX/2,
                                 worldSizeY/2,
                                 TargetSizeZ/2);   //its size

  G4LogicalVolume* logicTarget =
      new G4LogicalVolume(targetSolid,  //its solid
                          fpWaterMaterial,  //its material
                          "Target");    //its name

  new G4PVPlacement(0,  //no rotation
                    G4ThreeVector(),  //at (0,0,0)
                    "Target",  //its name
                    logicTarget,  //its logical volume
                    physiWorld,  //its mother  volume
                    false,        //no boolean operation
                    0);           //copy number

  // Visualization attributes
  
  G4VisAttributes* worldVisAtt =
      new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(worldVisAtt);

  G4VisAttributes* worldVisAtt1 = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  worldVisAtt1->SetVisibility(true);
  logicTarget->SetVisAttributes(worldVisAtt1);

  // Create Target G4Region and add logical volume

  fpRegion = new G4Region("Target");

  G4ProductionCuts* cuts = new G4ProductionCuts();

  G4double defCut = 1*nanometer;
  cuts->SetProductionCut(defCut,"gamma");
  cuts->SetProductionCut(defCut,"e-");
  cuts->SetProductionCut(defCut,"e+");
  cuts->SetProductionCut(defCut,"proton");

  fpRegion->SetProductionCuts(cuts);
  fpRegion->AddRootLogicalVolume(logicTarget); 

  return physiWorld;
}
