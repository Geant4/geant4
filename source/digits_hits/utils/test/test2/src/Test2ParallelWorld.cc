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
// $Id: Test2ParallelWorld.cc,v 1.3 2010-09-02 11:10:30 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Test2ParallelWorld.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"


Test2ParallelWorld::Test2ParallelWorld(G4String worldName)
  :G4VUserParallelWorld(worldName), fbConstructed(false)
{;}

Test2ParallelWorld::~Test2ParallelWorld()
{;}

void Test2ParallelWorld::Construct() {

  if(!fbConstructed)
  { 
    fbConstructed = true;
    SetupGeometry();
    SetupDetectors();
  }
}

void Test2ParallelWorld::SetupGeometry()
{
  //     
  // World
  //
  G4VPhysicalVolume * ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

  //                               
  // 3D nested phantom
  //
  // parameters
  G4double phantomSize[3] = {1.*m, 1.*m, 1.*m};
  G4int nSegment[3] = {30, 30, 30};

  //  phantom box
  G4VSolid * phantomSolid = new G4Box("PhantomBoxParallel", phantomSize[0], 
				      phantomSize[1], phantomSize[2]);
  G4LogicalVolume* phantomLogical = new G4LogicalVolume(phantomSolid, 0,
							"Phantom");
  new G4PVPlacement(0, G4ThreeVector(), phantomLogical, "PhantomParallel",
		    worldLogical, false, 0);

  G4String layerName[3] = {"layerX", "layerY", "layerZ"};
  G4VSolid * layerSolid[3];


  // replication along X
  layerSolid[0] = new G4Box(layerName[0],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1],
			    phantomSize[2]);
  fLayerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  new G4PVReplica(layerName[0], fLayerLogical[0], phantomLogical, kXAxis,
		  nSegment[0], phantomSize[0]/nSegment[0]*2.);

  // replication along Y
  layerSolid[1] = new G4Box(layerName[1],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]);
  fLayerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  new G4PVReplica(layerName[1], fLayerLogical[1], fLayerLogical[0], kYAxis,
		  nSegment[1], phantomSize[1]/nSegment[1]*2.);

  // replication along Z
  layerSolid[2] = new G4Box(layerName[2],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]/nSegment[2]);
  fLayerLogical[2] = new G4LogicalVolume(layerSolid[2], 0, layerName[2]);
  fPhantomPhys = new G4PVDivision(layerName[2], fLayerLogical[2],
				  fLayerLogical[1], kZAxis,
				  nSegment[2], 0.);

  //                                        
  // Visualization attributes
  //
  //  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes * invisatt = new G4VisAttributes(G4Colour(.0,.0,.0));
  invisatt->SetVisibility(false);
  phantomLogical->SetVisAttributes(invisatt);
  fLayerLogical[0]->SetVisAttributes(invisatt);
  fLayerLogical[1]->SetVisAttributes(invisatt);
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(1.,.7,.7));
  visatt->SetVisibility(true);
  fLayerLogical[2]->SetVisAttributes(visatt);

}


#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4SDParticleFilter.hh"

void Test2ParallelWorld::SetupDetectors() {
  //
  // primitive scorers
  //
  G4SDManager * sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(1);

  G4MultiFunctionalDetector* det = new G4MultiFunctionalDetector("ScoringWorld");

  // filters
  G4String filterName, particleName;
  G4SDParticleFilter* gammaFilter
    = new G4SDParticleFilter(filterName="gammaFilter",particleName="gamma");
  G4SDParticleFilter* electronFilter
    = new G4SDParticleFilter(filterName="electronFilter",particleName="e-");
  G4SDParticleFilter* positronFilter
    = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");

  // primitive scorers
  G4VPrimitiveScorer* primitive;
  primitive = new G4PSEnergyDeposit("eDep");
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthGamma");
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthElec");
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSTrackLength("trackLengthPosi");
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepGamma");
  primitive->SetFilter(gammaFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepElec");
  primitive->SetFilter(electronFilter);
  det->RegisterPrimitive(primitive);
  primitive = new G4PSNofStep("nStepPosi");
  primitive->SetFilter(positronFilter);
  det->RegisterPrimitive(primitive);

  sdManager->AddNewDetector(det);
  fLayerLogical[2]->SetSensitiveDetector(det);
  sdManager->SetVerboseLevel(0);
}

