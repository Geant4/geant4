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
// $Id$
//
// 

#include "Test2GeometryConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

Test2GeometryConstruction::Test2GeometryConstruction(G4double boxsize[3],
						     G4Material* material,
						     G4int segment[3])
  :fMaterial(material)
{
  phantomSize[0]=boxsize[0];
  phantomSize[1]=boxsize[1];
  phantomSize[2]=boxsize[2];

  nSegment[0]=segment[0];
  nSegment[1]=segment[1];
  nSegment[2]=segment[2];
}

Test2GeometryConstruction::~Test2GeometryConstruction()
{;}

G4VPhysicalVolume* 
Test2GeometryConstruction::ConstructGeometry(G4VPhysicalVolume* worldPhys)
{
  G4LogicalVolume* worldLogical = worldPhys->GetLogicalVolume();

  //  phantom box
  G4VSolid* phantomSolid = new G4Box("PhantomBox", phantomSize[0], phantomSize[1], phantomSize[2]);
  G4LogicalVolume* phantomLogical = new G4LogicalVolume(phantomSolid, fMaterial, "Phantom");
  G4cout << "logical" <<G4endl;
  new G4PVPlacement(0,G4ThreeVector(), phantomLogical, "Phantom",
                         worldLogical, false, 0);

  G4cout << phantomSize[0]/mm <<","<<phantomSize[1]/mm<<"," <<phantomSize[2]/mm <<G4endl;
  G4cout << nSegment[0]<<","<<nSegment[1] <<","<<nSegment[2] <<G4endl;

  G4String layerName[3] = {"layerX", "layerY", "layerZ"};
  G4VSolid * layerSolid[3];


  // replication along X
  layerSolid[0] = new G4Box(layerName[0],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1],
			    phantomSize[2]);
  fLayerLogical[0] = new G4LogicalVolume(layerSolid[0], fMaterial, layerName[0]);
  new G4PVReplica(layerName[0], fLayerLogical[0], phantomLogical, kXAxis,
		  nSegment[0], phantomSize[0]/nSegment[0]*2.);

  // replication along Y
  layerSolid[1] = new G4Box(layerName[1],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]);
  fLayerLogical[1] = new G4LogicalVolume(layerSolid[1], fMaterial, layerName[1]);
  new G4PVReplica(layerName[1], fLayerLogical[1], fLayerLogical[0], kYAxis,
		  nSegment[1], phantomSize[1]/nSegment[1]*2.);

  // replication along Z
  layerSolid[2] = new G4Box(layerName[2],
			    phantomSize[0]/nSegment[0],
			    phantomSize[1]/nSegment[1],
			    phantomSize[2]/nSegment[2]);
  fLayerLogical[2] = new G4LogicalVolume(layerSolid[2], fMaterial, layerName[2]);
  fPhantomPhys = new G4PVReplica(layerName[2], fLayerLogical[2], fLayerLogical[1], kZAxis,
  		  nSegment[2], phantomSize[2]/nSegment[2]*2.);
  //fPhantomPhys = new G4PVDivision(layerName[2], fLayerLogical[2], fLayerLogical[1], kZAxis,
  //		  nSegment[2], 0.);

  //                                        
  // Visualization attributes
  //
  //  worldLogical->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  phantomLogical->SetVisAttributes(simpleBoxVisAtt);

  G4VisAttributes * invisatt = new G4VisAttributes(G4Colour(.0,.0,.0));
  invisatt->SetVisibility(false);
  fLayerLogical[0]->SetVisAttributes(invisatt);
  fLayerLogical[1]->SetVisAttributes(invisatt);
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.8,.8,.8));
  visatt->SetVisibility(true);
  fLayerLogical[2]->SetVisAttributes(visatt);

  return worldPhys;

}


