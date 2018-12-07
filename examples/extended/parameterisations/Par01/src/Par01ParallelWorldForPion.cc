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
/// \file Par01/src/Par01ParallelWorldForPion.cc
/// \brief Implementation of the Par01ParallelWorldForPion class
//
//
#include "Par01ParallelWorldForPion.hh"
#include "Par01PionShowerModel.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"

#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01ParallelWorldForPion::Par01ParallelWorldForPion(G4String worldName)
  : G4VUserParallelWorld(worldName)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par01ParallelWorldForPion::~Par01ParallelWorldForPion()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01ParallelWorldForPion::Construct()
{
  // -- Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  
  // Build materials
  G4Material* air = nistManager->FindOrBuildMaterial("G4_AIR");
  
  // -------------------------------
  //  Build parallel/ghost geometry:
  // -------------------------------
  
  // -- Obtain clone of mass geometry world from GetWorld() base class utility:
  G4VPhysicalVolume* ghostWorld = GetWorld();
  // -- Why needed ? No default ?
  ghostWorld->GetLogicalVolume()->SetMaterial(air);
  G4Region* parallelWorldRegion = new G4Region("ParallelWorldRegion");
  parallelWorldRegion->AddRootLogicalVolume(ghostWorld->GetLogicalVolume());
  
  // -- Create box for encompassing Elec.+Had. calorimeters together:
  G4double detectSize = 125*cm;
  G4Box *ghostBox
    = new G4Box("GhostBox", detectSize+5*cm, detectSize+5*cm, 80*cm);
  // -- Build the subsequent logical volume:
  G4LogicalVolume* ghostLogical
    = new G4LogicalVolume(ghostBox,
                          air,               // -- no material : IMPOSSIBLE !!!
                          "GhostLogical", 
                          0, 0, 0);
  // -- And place this logical volume in the parallel geometry:
  new G4PVPlacement(0,G4ThreeVector(0., 0., 175*cm),
                    "GhostPhysical",
                    ghostLogical,
                    ghostWorld,false,0);
  
  // -----------------------
  //  Setup fast simulation:
  // -----------------------
  
  // -- Make Elec+Had logical volume becoming a region:
  G4cout << "Declaring region..." << G4endl;
  G4Region* ghostRegion = new G4Region("GhostCalorimeterRegion");
  G4cout << "... region declared. Setting it the ghostLogical volume..." << G4endl;
  ghostRegion->AddRootLogicalVolume(ghostLogical);
  G4cout << "... succes.\n" << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par01ParallelWorldForPion::ConstructSD()
{
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  
  G4Region* ghostRegion = regionStore->GetRegion("GhostCalorimeterRegion");
  // -- Attach fast simulation model (create the G4FastSimulationManager if needed):
  new Par01PionShowerModel("ghostPionShowerModel",ghostRegion);
  
}
