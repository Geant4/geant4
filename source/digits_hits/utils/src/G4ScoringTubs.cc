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
// $Id: G4ScoringTubs.cc,v 1.2 2007-08-22 01:40:29 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScoringTubs.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"
//#include "G4ScoringTubsParameterisation.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"


G4ScoringTubs::G4ScoringTubs(G4String wName)
  :G4VScoringMesh(wName), fSegmentDirection(-1), fRotationMatrix(0),
   fMeshElementLogical(0)
{
  fShape = cylinderMesh;

  fSize[0] = 0.*cm;
  fSize[1] = 1.*cm;
  fSize[2] = 1.*cm;
  fCenterPosition[0] = fCenterPosition[1] = fCenterPosition[2] = 0.*cm;
  fNSegment[0] = fNSegment[2] = 1;
  fNSegment[1] = 10;
}

G4ScoringTubs::~G4ScoringTubs()
{
}

void G4ScoringTubs::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {

    G4cerr << fWorldPhys->GetName() << G4endl;
    G4Exception(fWorldName+" has already been built.");

  } else {
    fConstructed = true;

    SetupGeometry(fWorldPhys);
  }
}



void G4ScoringTubs::SetupGeometry(G4VPhysicalVolume * fWorldPhys) {

  // World
  G4VPhysicalVolume * scoringWorld = fWorldPhys;
  G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  G4String tubsName("tubsMesh");
  if(fScoringMeshName.size() > 0) tubsName = fScoringMeshName;

  G4VSolid * tubsSolid = new G4Tubs(tubsName, fSize[0], fSize[1], fSize[2],
				    0., 360.*deg);
  G4LogicalVolume *  tubsLogical = new G4LogicalVolume(tubsSolid, 0, tubsName);
  new G4PVPlacement(fRotationMatrix, G4ThreeVector(fCenterPosition[0],
						   fCenterPosition[1],
						   fCenterPosition[2]),
		    tubsLogical, tubsName, worldLogical, false, 0);


  G4double fsegParam[3][3];
  G4int segOrder[3];
  if(fSegmentPositions.size() == 0) fSegmentDirection = 1; // r
  GetSegmentOrder(fSegmentDirection, fNSegment, segOrder, fsegParam);
  EAxis axis[3] = {kRho, kPhi, kZAxis};



  G4String layerName[2] = {tubsName + "_nest1",  tubsName + "_nest2"};
  G4VSolid * layerSolid[2]; 
  G4LogicalVolume * layerLogical[2];

  // fisrt nested layer
  //G4cout << "layer 1 :" << G4endl;
  layerSolid[0] = new G4Tubs(layerName[0],
			     fSize[0],
			     fSize[0]+fsegParam[0][0],
			     fsegParam[0][2],
			     0., fsegParam[0][1]);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  if(fNSegment[segOrder[0]] > 1) {
    G4double width = fsegParam[0][segOrder[0]];
    if(axis[segOrder[0]] == kZAxis) width *= 2.;
    new G4PVReplica(layerName[0], layerLogical[0], tubsLogical, axis[segOrder[0]],
		    fNSegment[segOrder[0]], width, 0.);
  } else if(fNSegment[segOrder[0]] == 1) {
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[0], layerName[0], tubsLogical, false, 0);
  } else {
    G4cerr << "G4ScoringTubs::SetupGeometry() : invalid parameter ("
	   << fNSegment[segOrder[0]] << "," << segOrder[0] << ") "
	   << "in placement of the first nested layer." << G4endl;
  }

  /*
  G4cout << fSize[0] << ", " << fsegParam[0][0] << ", "
	 << fsegParam[0][2] << ", " << fsegParam[0][1]/deg << G4endl;
  G4cout << fNSegment[segOrder[0]] << ", " <<  fsegParam[0][segOrder[0]] << G4endl;
  */

  // second nested layer
  //G4cout << "layer 2 :" << G4endl;
  layerSolid[1] = new G4Tubs(layerName[1],
			     fSize[0],
			     fSize[0]+fsegParam[1][0],
			     fsegParam[1][2],
			     0., fsegParam[1][1]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  if(fNSegment[segOrder[1]] > 1)  {
    G4double width = fsegParam[1][segOrder[1]];
    if(axis[segOrder[1]] == kZAxis) width *= 2.;
    new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], axis[segOrder[1]],
		    fNSegment[segOrder[1]], width, 0.);
  } else if(fNSegment[segOrder[1]] == 1)
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[1], layerName[1], layerLogical[0], false, 0);
  else
    G4cerr << "G4ScoringTubs::SetupGeometry() : invalid parameter ("
	   << fNSegment[segOrder[1]] << "," << segOrder[1] << ") "
	   << "in placement of the second nested layer." << G4endl;


  /*
  G4cout << fSize[0] << ", " << fsegParam[1][0] << ", "
	 << fsegParam[1][2] << ", " << fsegParam[1][1]/deg << G4endl;
  G4cout << fNSegment[segOrder[1]] << ", " <<  fsegParam[1][segOrder[1]] << G4endl;
  */

  // mesh elements
  //G4cout << "mesh elements :" << G4endl;
  G4String elementName = tubsName +"_element";
  G4VSolid * elementSolid = new G4Tubs(elementName,
				       fSize[0],
				       fSize[0]+fsegParam[2][0],
				       fsegParam[2][2],
				       0., fsegParam[2][1]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, 0, elementName);
  if(fNSegment[segOrder[2]] > 1) {
    /*
    if(fSegmentPositions.size() > 0) {
      G4double motherDims[3] ={fSize[0]/fsegParam[2][0],
			       fSize[1]/fsegParam[2][1],
			       fSize[2]/fsegParam[2][2]};
      G4int nelement = fSegmentPositions.size() + 1;
      //G4ScoringTubsParameterisation * param =
      G4VPVParameterisation * param =
	new G4ScoringTubsParameterisation(axis[2], motherDims, fSegmentPositions);
      new G4PVParameterised(elementName,
			    fMeshElementLogical,
			    layerLogical[1],
			    axis[2],
			    nelement,
			    param);
    } else {
    */
    G4double width = fsegParam[2][segOrder[2]];
    if(axis[segOrder[2]] == kZAxis) width *= 2.;
    new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], axis[segOrder[2]],
		    fNSegment[segOrder[2]], width, 0.);
    //}
  } else if(fNSegment[segOrder[2]] == 1) {
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), fMeshElementLogical, elementName, layerLogical[1], false, 0);
  } else {
    G4cerr << "G4ScoringTubs::SetupGeometry() : "
	   << "invalid parameter (" << fNSegment[segOrder[2]] << ") "
	   << "in mesh element placement." << G4endl;
  }

  /*
  G4cout << fSize[0] << ", " << fsegParam[2][0] << ", "
	 << fsegParam[2][2] << ", " << fsegParam[2][1]/deg << G4endl;
  G4cout << fNSegment[segOrder[2]] << ", " <<  fsegParam[2][segOrder[2]] << G4endl;
  */

  //
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5));
  visatt->SetVisibility(true);
  fMeshElementLogical->SetVisAttributes(visatt);
}

void G4ScoringTubs::RegisterPrimitives(std::vector<G4VPrimitiveScorer *> & vps) {


  if(fMFD == 0) fMFD = new G4MultiFunctionalDetector(fWorldName);

  std::vector<G4VPrimitiveScorer *>::iterator itr = vps.begin();
  for(; itr != vps.end(); itr++) {
    fMFD->RegisterPrimitive(*itr);

    fPS.push_back(*itr);
  }

  G4SDManager::GetSDMpointer()->AddNewDetector(fMFD);
  fMeshElementLogical->SetSensitiveDetector(fMFD);
}

void G4ScoringTubs::List() const {
  G4cout << "G4ScoringTubs" << G4endl;
}


void G4ScoringTubs::GetSegmentOrder(G4int segDir, G4int nseg[3], 
				    G4int segOrd[3], G4double segParam[3][3]) {

  segParam[2][0] = (fSize[1] - fSize[0])/nseg[0]; // r
  segParam[2][1] = 360.*deg/nseg[1]; // phi
  segParam[2][2] = fSize[2]/nseg[2]; // z

  if(segDir == -1 || segDir == 1) { // in r direction, it segments z -> phi -> r directions by turns.
    segOrd[0] = 2;
    segOrd[1] = 1;
    segOrd[2] = 0;
    segParam[0][0] = fSize[1] - fSize[0];
    segParam[0][1] = 360.*deg;
    segParam[0][2] = fSize[2]/nseg[2];
    segParam[1][0] = fSize[1] - fSize[0];
    segParam[1][1] = 360.*deg/nseg[1];
    segParam[1][2] = fSize[2]/nseg[2];

  } else if(segDir == 3) { // in z direction, it segments phi -> r -> z directions by turns.
    segOrd[0] = 1;
    segOrd[1] = 0;
    segOrd[2] = 2;
    segParam[0][0] = fSize[1] - fSize[0];
    segParam[0][1] = 360.*deg/nseg[1];
    segParam[0][2] = fSize[2];
    segParam[1][0] = (fSize[1] - fSize[0])/nseg[0];
    segParam[1][1] = 360.*deg/nseg[1];
    segParam[1][2] = fSize[2];

  } else if(segDir == 2) { // in phi direction, it segments z -> r -> phi directions by turns.
    segOrd[0] = 2;
    segOrd[1] = 0;
    segOrd[2] = 1;
    segParam[0][0] = fSize[1] - fSize[0];
    segParam[0][1] = 360.*deg;
    segParam[0][2] = fSize[2]/nseg[2];
    segParam[1][0] = (fSize[1] - fSize[0])/nseg[0];
    segParam[1][1] = 360.*deg;
    segParam[1][2] = fSize[2]/nseg[2];

  } else {
    G4cerr << "G4ScoringTubs : " 
	   << " The segment direction (" << segDir
	   << ") is invalid value." << G4endl;
    segOrd[0] = 0;
    segOrd[1] = 1;
    segOrd[2] = 2;
    segParam[0][0] = segParam[0][1] = segParam[0][2] = 
      segParam[1][0] = segParam[1][1] = segParam[1][2] = 
      segParam[2][0] = segParam[2][1] = segParam[2][2] =  1.;
  }

}

