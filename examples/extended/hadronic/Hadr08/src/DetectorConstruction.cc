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
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4RegionStore.hh"
#include "BiasingOperator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
  fSiliconTrackerLogical( nullptr ), fEmCaloLogical( nullptr ), fHadCaloLogical( nullptr ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {

  G4cout << "\nDetectorConstruction....\n" << G4endl;
  
  //--------- Material definition ---------
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* air     = nistManager->FindOrBuildMaterial( "G4_AIR" );
  G4Material* csi     = nistManager->FindOrBuildMaterial( "G4_CESIUM_IODIDE" );
  G4Material* silicon = nistManager->FindOrBuildMaterial( "G4_Si" );
  G4Material* iron    = nistManager->FindOrBuildMaterial( "G4_Fe" );

  // World
  G4Box* worldBox = new G4Box( "WorldBox", 200.0*cm, 200.0*cm, 200.0*cm );
  G4LogicalVolume* worldLogical = new G4LogicalVolume( worldBox, air, "WorldLogical", 0, 0, 0 );
  G4PVPlacement* worldPhys = new G4PVPlacement( 0, G4ThreeVector(), "WorldPhysical", 
                                                worldLogical, 0, false, 0 );

  // Silicon Tracker: for simplicity and to have more hadronic interactions,
  //                  it is single layer of silicon (20 cm thick)
  G4double halfDetectorXYsize = 100.0*cm;
  G4Box* siliconBox = new G4Box( "siliconBox", halfDetectorXYsize, halfDetectorXYsize, 10.0*cm );
  fSiliconTrackerLogical = new G4LogicalVolume( siliconBox, silicon, 
                                                "SiliconTrackerLogical", 0, 0, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.0 ), "SiliconTrackerPhysical", 
                    fSiliconTrackerLogical, worldPhys, false, 0 );
  
  // Electromagnetic (EM) calorimeter: an homogeneous crystal (40 cm thick)
  G4Box* emCaloBox = new G4Box( "EmCaloBox", halfDetectorXYsize, halfDetectorXYsize, 20.0*cm );
  fEmCaloLogical = new G4LogicalVolume( emCaloBox, csi, "EmCalorimeterLogical", 0, 0, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 35.0*cm ), "EmCalorimeterPhysical",
                     fEmCaloLogical, worldPhys, false, 0 );

  // The hadronic (HAD) calorimeter: an homogeneous block of metal (100 cm thick)
  G4Box* hadCaloBox = new G4Box( "HadCaloBox", halfDetectorXYsize, halfDetectorXYsize, 50.0*cm );
  fHadCaloLogical = new G4LogicalVolume( hadCaloBox, iron, "HadCaloLogical", 0, 0, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 110.0*cm ), "HadCaloPhysical", 
                     fHadCaloLogical, worldPhys, false, 0 );

  // Define the tracker region
  G4Region* trackerRegion = new G4Region( "Tracker_region" );
  trackerRegion->AddRootLogicalVolume( fSiliconTrackerLogical );
  std::vector< double > cuts;
  cuts.push_back( 10.0*cm ); cuts.push_back( 10.0*cm ); 
  cuts.push_back( 10.0*cm ); cuts.push_back( 10.0*cm );
  trackerRegion->SetProductionCuts( new G4ProductionCuts );
  trackerRegion->GetProductionCuts()->SetProductionCuts( cuts );

  // Define the EM calorimeter region
  G4Region* emCaloRegion = new G4Region( "EM_calo_region" );
  emCaloRegion->AddRootLogicalVolume( fEmCaloLogical );
  cuts.clear(); 
  cuts.push_back( 20.0*cm ); cuts.push_back( 20.0*cm );
  cuts.push_back( 20.0*cm ); cuts.push_back( 20.0*cm );
  emCaloRegion->SetProductionCuts( new G4ProductionCuts );
  emCaloRegion->GetProductionCuts()->SetProductionCuts( cuts );

  // Define the HAD calorimeter region
  G4Region* hadCaloRegion = new G4Region( "HAD_calo_region" );
  hadCaloRegion->AddRootLogicalVolume( fHadCaloLogical );
  cuts.clear();
  cuts.push_back( 50.0*cm ); cuts.push_back( 50.0*cm ); 
  cuts.push_back( 50.0*cm ); cuts.push_back( 50.0*cm );
  hadCaloRegion->SetProductionCuts( new G4ProductionCuts );
  hadCaloRegion->GetProductionCuts()->SetProductionCuts( cuts );

  return worldPhys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {
  // Instanciate the biasing operator and specify the particles of interest (i.e. p, n, pi+, pi-)
  // and the logical volume(s) where this operator should be applied. 
  BiasingOperator* biasingOperator = new BiasingOperator;
  biasingOperator->AddParticle( "proton" );
  biasingOperator->AddParticle( "neutron" );
  biasingOperator->AddParticle( "pi+" );
  biasingOperator->AddParticle( "pi-" );
  biasingOperator->AttachTo( fSiliconTrackerLogical );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
