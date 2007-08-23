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
// $Id: G4ScoringBox.cc,v 1.9 2007-08-23 02:21:22 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScoringBox.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"
#include "G4ScoringBoxParameterisation.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"


G4ScoringBox::G4ScoringBox(G4String wName)
  :G4VScoringMesh(wName), fRotationMatrix(0), fSegmentDirection(-1),
   fMeshElementLogical(0)
{
  fShape = boxMesh;

  fSize[0] = fSize[1] = fSize[2] = 1.*cm;
  fCenterPosition[0] = fCenterPosition[1] = fCenterPosition[2] = 0.*cm;
  fNSegment[0] = fNSegment[1] = fNSegment[2] = 1;
}

G4ScoringBox::~G4ScoringBox()
{
}

void G4ScoringBox::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {

    G4cerr << fWorldPhys->GetName() << G4endl;
//////////////////////////////////////    G4Exception(fWorldName+" has already been built.");

  } else {
    fConstructed = true;

    SetupGeometry(fWorldPhys);
  }
}



void G4ScoringBox::SetupGeometry(G4VPhysicalVolume * fWorldPhys) {

  //G4cout << "G4ScoringBox::SetupGeometry() ..." << G4endl;

  // World
  G4VPhysicalVolume * scoringWorld = fWorldPhys;
  G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  //G4cout << "box mesh" << G4endl;
  G4String boxName("boxMesh");
  if(fScoringMeshName.size() > 0) boxName = fScoringMeshName;

  G4VSolid * boxSolid = new G4Box(boxName, fSize[0], fSize[1], fSize[2]);
  G4LogicalVolume *  boxLogical = new G4LogicalVolume(boxSolid, 0, boxName);
  new G4PVPlacement(fRotationMatrix, G4ThreeVector(fCenterPosition[0],
						   fCenterPosition[1],
						   fCenterPosition[2]),
		    boxLogical, boxName, worldLogical, false, 0);
  //G4cout << fSize[0] << ", " << fSize[1] << ", " << fSize[2] << G4endl;

  G4double fsegment[3][3];
  G4int segOrder[3];
  GetSegmentOrder(fSegmentDirection, fNSegment, segOrder, fsegment);
  EAxis axis[3] = {kXAxis, kYAxis, kZAxis};



  G4String layerName[2] = {boxName + "_nest1",  boxName + "_nest2"};
  G4VSolid * layerSolid[2]; 
  G4LogicalVolume * layerLogical[2];

  // fisrt nested layer
  //G4cout << "layer 1 :" << G4endl;
  layerSolid[0] = new G4Box(layerName[0],
			    fSize[0]/fsegment[0][0],
			    fSize[1]/fsegment[0][1],
			    fSize[2]/fsegment[0][2]);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  if(fNSegment[segOrder[0]] > 1) 
    new G4PVReplica(layerName[0], layerLogical[0], boxLogical, axis[segOrder[0]],
		    fNSegment[segOrder[0]], fSize[segOrder[0]]/fsegment[0][segOrder[0]]*2.);
  else if(fNSegment[segOrder[0]] == 1)
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[0], layerName[0], boxLogical, false, 0);
  else
    G4cerr << "G4ScoringBox::SetupGeometry() : invalid parameter ("
	   << fNSegment[segOrder[0]] << "," << segOrder[0] << ") "
	   << "in placement of the first nested layer." << G4endl;

  /*
  G4cout << fSize[0]/fsegment[0][0] << ", "
	 << fSize[1]/fsegment[0][1] << ", "
	 << fSize[2]/fsegment[0][2] << G4endl;
  G4cout << layerName[0] << ": " << axis[segOrder[0]] << ", "
	 << fNSegment[segOrder[0]] << ", "
	 << fSize[segOrder[0]]/fsegment[0][segOrder[0]] << G4endl;
  */

  // second nested layer
  //G4cout << "layer 2 :" << G4endl;
  layerSolid[1] = new G4Box(layerName[1],
			    fSize[0]/fsegment[1][0],
			    fSize[1]/fsegment[1][1],
			    fSize[2]/fsegment[1][2]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  if(fNSegment[segOrder[1]] > 1) 
    new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], axis[segOrder[1]],
		    fNSegment[segOrder[1]], fSize[segOrder[1]]/fsegment[1][segOrder[1]]*2.);
  else if(fNSegment[segOrder[1]] == 1)
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[1], layerName[1], layerLogical[0], false, 0);
  else
    G4cerr << "G4ScoringBox::SetupGeometry() : invalid parameter ("
	   << fNSegment[segOrder[1]] << "," << segOrder[1] << ") "
	   << "in placement of the second nested layer." << G4endl;

  /*
  G4cout << fSize[0]/fsegment[1][0] << ", "
	 << fSize[1]/fsegment[1][1] << ", "
	 << fSize[2]/fsegment[1][2] << G4endl;
  G4cout << layerName[1] << ": " << axis[segOrder[1]] << ", "
	 << fNSegment[segOrder[1]] << ", "
	 << fSize[segOrder[1]]/fsegment[1][segOrder[1]] << G4endl;
  */

  // mesh elements
  //G4cout << "mesh elements :" << G4endl;
  G4String elementName = boxName +"_element";
  G4VSolid * elementSolid = new G4Box(elementName, fSize[0]/fsegment[2][0],
				      fSize[1]/fsegment[2][1], fSize[2]/fsegment[2][2]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, 0, elementName);
  if(fNSegment[segOrder[2]] > 1) 
    if(fSegmentPositions.size() > 0) {
      G4double motherDims[3] ={fSize[0]/fsegment[1][0],
			       fSize[1]/fsegment[1][1],
			       fSize[2]/fsegment[1][2]};
      G4int nelement = fNSegment[segOrder[2]];
      fSegmentPositions.push_back(fSize[segOrder[2]]*2.);
      //G4ScoringBoxParameterisation * param =
      G4VPVParameterisation * param =
	new G4ScoringBoxParameterisation(axis[segOrder[2]], motherDims, fSegmentPositions);
      new G4PVParameterised(elementName,
			    fMeshElementLogical,
			    layerLogical[1],
			    axis[segOrder[2]],
			    nelement,
			    param);

      /*
      G4cout << motherDims[0] << ", " << motherDims[1] << ", " << motherDims[2] << G4endl;
      for(int i = 0; i < (int)fSegmentPositions.size(); i++)
	G4cout << fSegmentPositions[i] << ", ";
      G4cout << G4endl;
      */

    } else {
      new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], axis[segOrder[2]],
		      fNSegment[segOrder[2]], fSize[segOrder[2]]/fsegment[2][segOrder[2]]*2.);
    }
  else if(fNSegment[segOrder[2]] == 1)
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), fMeshElementLogical, elementName, layerLogical[1], false, 0);
  else
    G4cerr << "G4ScoringBox::SetupGeometry() : "
	   << "invalid parameter (" << fNSegment[segOrder[2]] << ") "
	   << "in mesh element placement." << G4endl;

  /*
  G4cout << fSize[0]/fsegment[2][0] << ", "
	 << fSize[1]/fsegment[2][1] << ", "
	 << fSize[2]/fsegment[2][2] << G4endl;
  G4cout << elementName << ": " << axis[segOrder[2]] << ", "
	 << fNSegment[segOrder[2]] << ", "
	 << fSize[segOrder[2]]/fsegment[2][segOrder[2]] << G4endl;
  */

  //
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5));
  visatt->SetVisibility(true);
  fMeshElementLogical->SetVisAttributes(visatt);
}

void G4ScoringBox::RegisterPrimitives(std::vector<G4VPrimitiveScorer *> & vps) {


  if(fMFD == 0) fMFD = new G4MultiFunctionalDetector(fWorldName);

  std::vector<G4VPrimitiveScorer *>::iterator itr = vps.begin();
  for(; itr != vps.end(); itr++) {
    fMFD->RegisterPrimitive(*itr);

    fPS.push_back(*itr);
  }

  G4SDManager::GetSDMpointer()->AddNewDetector(fMFD);
  fMeshElementLogical->SetSensitiveDetector(fMFD);
}

void G4ScoringBox::List() const {
  G4cout << "G4ScoringBox" << G4endl;
}


void G4ScoringBox::GetSegmentOrder(G4int segDir, G4int nseg[3], G4int segOrd[3], G4double segfact[3][3]) {

  if(segDir == -1 || segDir == 1) { // in x direction, it segments z -> y -> x directions by turns.
    segOrd[0] = 2;
    segOrd[1] = 1;
    segOrd[2] = 0;
    segfact[0][0] = 1.;
    segfact[0][1] = 1.;
    segfact[0][2] = G4double(nseg[2]);
    segfact[1][0] = 1.;
    segfact[1][1] = G4double(nseg[1]);
    segfact[1][2] = G4double(nseg[2]);
    segfact[2][0] = G4double(nseg[0]);
    segfact[2][1] = G4double(nseg[1]);
    segfact[2][2] = G4double(nseg[2]);

  } else if(segDir == 2) { // in y direction, it segments x -> z -> y directions by turns.
    segOrd[0] = 0;
    segOrd[1] = 2;
    segOrd[2] = 1;
    segfact[0][0] = G4double(nseg[0]);
    segfact[0][1] = 1.;
    segfact[0][2] = 1.;
    segfact[1][0] = G4double(nseg[0]);
    segfact[1][1] = 1.;
    segfact[1][2] = G4double(nseg[2]);
    segfact[2][0] = G4double(nseg[0]);
    segfact[2][1] = G4double(nseg[1]);
    segfact[2][2] = G4double(nseg[2]);

  } else if(segDir == 3) { // in z direction, it segments y -> x -> z directions by turns.
    segOrd[0] = 1;
    segOrd[1] = 0;
    segOrd[2] = 2;
    segfact[0][0] = 1.;
    segfact[0][1] = G4double(nseg[1]);
    segfact[0][2] = 1.;
    segfact[1][0] = G4double(nseg[0]);
    segfact[1][1] = G4double(nseg[1]);
    segfact[1][2] = 1.;
    segfact[2][0] = G4double(nseg[0]);
    segfact[2][1] = G4double(nseg[1]);
    segfact[2][2] = G4double(nseg[2]);
  } else {
    G4cerr << "G4ScoringBox : " 
	   << " The segment direction (" << segDir
	   << ") is invalid value." << G4endl;
    segOrd[0] = 0;
    segOrd[1] = 1;
    segOrd[2] = 2;
    segfact[0][0] = segfact[0][1] = segfact[0][2] = 
      segfact[1][0] = segfact[1][1] = segfact[1][2] = 
      segfact[2][0] = segfact[2][1] = segfact[2][2] =  1.;
  }
}

