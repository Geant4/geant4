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
// $Id: A01ParallelWorld.cc,v 1.3 2006-12-13 15:49:20 gunter Exp $
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

#include "G4VisAttributes.hh"

A01ParallelWorld::A01ParallelWorld(G4String worldName)
:G4VUserParallelWorld(worldName)
{;}

A01ParallelWorld::~A01ParallelWorld() 
{;}

void A01ParallelWorld::Construct() 
{
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

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


