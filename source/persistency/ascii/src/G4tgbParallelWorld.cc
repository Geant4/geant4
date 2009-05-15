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
#include "G4tgbParallelWorld.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"

#include "G4tgrVolumeMgr.hh"
#include "G4tgbVolumeMgr.hh"
#include "G4tgbVolume.hh"
#include "G4tgbDetectorBuilder.hh"


G4tgbParallelWorld::G4tgbParallelWorld(G4String worldName, G4int index)
  : G4VUserParallelWorld(worldName), theIndex(index)
{;}

G4tgbParallelWorld::~G4tgbParallelWorld()
{;}

#include "G4PhysicalVolumeStore.hh"

void G4tgbParallelWorld::Construct()
{
  // -------------------------------
  //  Build parallel/ghost geometry:
  // -------------------------------
  G4VPhysicalVolume* ghostWorld = GetWorld();

  //  G4cout << " G4tgbParallelWorld::Construct() ghostWorld " << ghostWorld->GetName() << " " << ghostWorld << G4endl;
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  G4tgbVolumeMgr* tgbVolmgr = G4tgbVolumeMgr::GetInstance();
  tgbVolmgr->RegisterMe(worldLogical);

  //  G4cout << " G4tgbParallelWorld::Construct() worldLogical " << worldLogical<< G4endl;
  
  // -- Obtain clone of mass geometry world from GetWorld() base class utility:
  tgbVolmgr->RegisterMe(ghostWorld);
  G4tgbDetectorBuilder* detectorBuilder = tgbVolmgr->GetDetectorBuilder();
  const G4tgrVolume* tgrVoltop = G4tgrVolumeMgr::GetInstance()->GetTopVolume(theIndex);  
  detectorBuilder->ConstructDetector( tgrVoltop, theIndex );

    //check that there are not two

    // No need of Construct in parallelworld 
    //- set list of logical volumes that are placed in the world volume
 
  // -- Why needed ? No default ?
  //-  ghostWorld->GetLogicalVolume()->SetMaterial(Air);
  //  G4Region* parallelWorldRegion = new G4Region("ParallelWorldRegion");
  //  parallelWorldRegion->AddRootLogicalVolume(ghostWorld->GetLogicalVolume());
  /*
    
  // place volumes in the parallel world here. For example ...
  //
  G4Box * ghostSolid = new G4Box("GhostdBox", 1.*cm, 1.*cm, 1.*cm);
  G4LogicalVolume * ghostLogical
        = new G4LogicalVolume(ghostSolid, 0, "GhostLogical", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(), ghostLogical,
                    "GhostPhysical", worldLogical, 0, 0);

  */

  // get placements in world 
  const G4PhysicalVolumeStore* lvs = G4PhysicalVolumeStore::GetInstance();
  std::vector<G4VPhysicalVolume*>::const_iterator lvcite;
  for( lvcite = lvs->begin(); lvcite != lvs->end(); lvcite++ ) {
    if( (*lvcite)->GetMotherLogical() == worldLogical )
      {
	//	G4cout << " placement in world " << (*lvcite)->GetName() << G4endl;
      }
  }
}
