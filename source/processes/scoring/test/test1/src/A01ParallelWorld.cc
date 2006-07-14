//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: A01ParallelWorld.cc,v 1.1 2006-07-14 14:43:27 asaim Exp $
// --------------------------------------------------------------
//

#include "A01ParallelWorld.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSNofStep.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofSecondary.hh"
#include "G4PSEnergyDeposit.hh"

A01ParallelWorld::A01ParallelWorld(G4String worldName)
:G4VUserParallelWorld(worldName)
{;}

A01ParallelWorld::~A01ParallelWorld() 
{;}

void A01ParallelWorld::Construct() 
{
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
 
  G4Box * ghostSolid = new G4Box("GhostdBox", 60.*cm, 60.*cm, 60.*cm);
  G4LogicalVolume * ghostLogical
        = new G4LogicalVolume(ghostSolid, 0, "GhostLogical", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(), ghostLogical, "GhostPhysical",
                            worldLogical, 0, 0);
  G4Box * layerSolid = new G4Box("GhostLayerBox", 60.*cm, 60.*cm, 2.*cm);
  G4LogicalVolume * layerLogical
        = new G4LogicalVolume(layerSolid, 0, "GhostLayerLogical", 0, 0, 0);
  new G4PVReplica("GhostLayerPhysical",layerLogical,ghostLogical,kZAxis,30,4.*cm);


  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SDname;
  G4MultiFunctionalDetector* aSD = new G4MultiFunctionalDetector(SDname="ParallelWorld");
  SDman->AddNewDetector(aSD);
  layerLogical->SetSensitiveDetector(aSD);

  aSD->RegisterPrimitive(new G4PSNofStep("NofStep"));
  aSD->RegisterPrimitive(new G4PSTrackLength("TrackLength"));
  aSD->RegisterPrimitive(new G4PSNofSecondary("NofSecondary"));
  aSD->RegisterPrimitive(new G4PSEnergyDeposit("EnergyDeposit"));
}


