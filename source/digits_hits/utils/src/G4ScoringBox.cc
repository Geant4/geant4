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
// $Id: G4ScoringBox.cc,v 1.26 2007-09-03 13:32:16 akimura Exp $
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
#include "G4VVisManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit3D.hh"


G4ScoringBox::G4ScoringBox(G4String wName)
  :G4VScoringMesh(wName), fSegmentDirection(-1),
   fMeshElementLogical(0)
{
  fShape = boxMesh;
}

G4ScoringBox::~G4ScoringBox()
{
}

void G4ScoringBox::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {

    G4cerr << fWorldPhys->GetName() << G4endl;
//////////////////////////////////////    G4Exception(fWorldName+" has already been built.");
    ResetScore();

  } else {
    fConstructed = true;

    SetupGeometry(fWorldPhys);
  }
}



void G4ScoringBox::SetupGeometry(G4VPhysicalVolume * fWorldPhys) {

  if(verboseLevel > 10) G4cout << "G4ScoringBox::SetupGeometry() ..." << G4endl;

  // World
  G4VPhysicalVolume * scoringWorld = fWorldPhys;
  G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  if(verboseLevel > 10) G4cout << fWorldName << G4endl;
  G4String boxName = fWorldName;

  if(verboseLevel > 10) G4cout << fSize[0] << ", " << fSize[1] << ", " << fSize[2] << G4endl;
  G4VSolid * boxSolid = new G4Box(boxName+"0", fSize[0], fSize[1], fSize[2]);
  G4LogicalVolume *  boxLogical = new G4LogicalVolume(boxSolid, 0, boxName);
  new G4PVPlacement(fRotationMatrix, fCenterPosition,
		    boxLogical, boxName+"0", worldLogical, false, 0);

  //G4double fsegment[3][3];
  //G4int segOrder[3];
  //GetSegmentOrder(fSegmentDirection, fNSegment, segOrder, fsegment);
  //EAxis axis[3] = {kXAxis, kYAxis, kZAxis};



  G4String layerName[2] = {boxName + "1",  boxName + "2"};
  G4VSolid * layerSolid[2]; 
  G4LogicalVolume * layerLogical[2];

  // fisrt nested layer (replicated to x direction)
  if(verboseLevel > 10) G4cout << "layer 1 :" << G4endl;
  layerSolid[0] = new G4Box(layerName[0],
			    fSize[0]/fNSegment[0],
			    fSize[1],
			    fSize[2]);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  if(fNSegment[0] > 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Replicate to x direction" << G4endl;
    new G4PVReplica(layerName[0], layerLogical[0], boxLogical, kXAxis,
		    fNSegment[0], fSize[0]/fNSegment[0]*2.);
  } else if(fNSegment[0] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[0], layerName[0], boxLogical, false, 0);
  } else
    G4cerr << "G4ScoringBox::SetupGeometry() : invalid parameter ("
	   << fNSegment[0] << ") "
	   << "in placement of the first nested layer." << G4endl;

  if(verboseLevel > 10) {
    G4cout << fSize[0]/fNSegment[0] << ", "
	   << fSize[1] << ", "
	   << fSize[2] << G4endl;
    G4cout << layerName[0] << ": kXAxis, "
	   << fNSegment[0] << ", "
	   << 2.*fSize[0]/fNSegment[0] << G4endl;
  }

  // second nested layer (replicated to y direction)
  if(verboseLevel > 10) G4cout << "layer 2 :" << G4endl;
  layerSolid[1] = new G4Box(layerName[1],
			    fSize[0]/fNSegment[0],
			    fSize[1]/fNSegment[1],
			    fSize[2]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  if(fNSegment[1] > 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Replicate to y direction" << G4endl;
    new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], kYAxis,
		    fNSegment[1], fSize[1]/fNSegment[1]*2.);
  } else if(fNSegment[1] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[1], layerName[1], layerLogical[0], false, 0);
  } else
    G4cerr << "G4ScoringBox::SetupGeometry() : invalid parameter ("
	   << fNSegment[1] << ") "
	   << "in placement of the second nested layer." << G4endl;

  if(verboseLevel > 10) {
    G4cout << fSize[0]/fNSegment[0] << ", "
	   << fSize[1]/fNSegment[1] << ", "
	   << fSize[2] << G4endl;
    G4cout << layerName[1] << ": kYAxis, "
	   << fNSegment[1] << ", "
	   << 2.*fSize[1]/fNSegment[1] << G4endl;
  }

  // mesh elements (replicated to z direction)
  if(verboseLevel > 10) G4cout << "mesh elements :" << G4endl;
  G4String elementName = boxName +"3";
  G4VSolid * elementSolid = new G4Box(elementName,
				      fSize[0]/fNSegment[0],
				      fSize[1]/fNSegment[1],
				      fSize[2]/fNSegment[2]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, 0, elementName);
  if(fNSegment[2] > 1) 
    if(fSegmentPositions.size() > 0) {
      if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Parameterise to z direction" << G4endl;
      G4double motherDims[3] ={fSize[0]/fNSegment[0],
			       fSize[1]/fNSegment[1],
			       fSize[2]/fNSegment[2]};
      G4int nelement = fNSegment[2];
      fSegmentPositions.push_back(fSize[2]*2.);
      //G4ScoringBoxParameterisation * param =
      G4VPVParameterisation * param =
	new G4ScoringBoxParameterisation(kZAxis, motherDims, fSegmentPositions);
      new G4PVParameterised(elementName,
			    fMeshElementLogical,
			    layerLogical[1],
			    kZAxis,
			    nelement,
			    param);

      if(verboseLevel > 10) {
	G4cout << motherDims[0] << ", " << motherDims[1] << ", " << motherDims[2] << G4endl;
	for(int i = 0; i < (int)fSegmentPositions.size(); i++)
	  G4cout << fSegmentPositions[i] << ", ";
	G4cout << G4endl;
      }

    } else {
      if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Replicate to z direction" << G4endl;

      new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], kZAxis,
		      fNSegment[2], 2.*fSize[2]/static_cast<G4double>(fNSegment[2]));
    }
  else if(fNSegment[2] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), fMeshElementLogical, elementName, layerLogical[1], false, 0);
  } else
    G4cerr << "G4ScoringBox::SetupGeometry() : "
	   << "invalid parameter (" << fNSegment[2] << ") "
	   << "in mesh element placement." << G4endl;

  if(verboseLevel > 10) {
    G4cout << fSize[0]/fNSegment[0] << ", "
	   << fSize[1]/fNSegment[1] << ", "
	   << fSize[2]/fNSegment[2] << G4endl;
    G4cout << elementName << ": kZAxis, "
	   << fNSegment[2] << ", "
	   << 2.*fSize[2]/fNSegment[2] << G4endl;
  }


  // set the sensitive detector
  fMeshElementLogical->SetSensitiveDetector(fMFD);
  

  // vis. attributes
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5));
  visatt->SetVisibility(true);
  layerLogical[0]->SetVisAttributes(visatt);
  layerLogical[1]->SetVisAttributes(visatt);
  fMeshElementLogical->SetVisAttributes(visatt);
}


void G4ScoringBox::List() const {
  G4cout << "G4ScoringBox : " << fWorldName << G4endl;
  G4cout << " Shape: Box mesh" << G4endl;

  G4VScoringMesh::List();

  G4cout << "# of G4THitsMap : " << fMap.size() << G4endl;

  std::map<G4String, G4THitsMap<G4double>* >::const_iterator itr = fMap.begin();
  for(; itr != fMap.end(); itr++) {
    G4cout << "[" << itr->first << "]" << G4endl;
    /*
    G4THitsMap<G4double> * map = itr->second;
    std::map<G4int, G4double*>::iterator itrMap = map->GetMap()->begin();
    G4int q[3], idx;
    G4double tot = 0.;
    G4cout << "    position       value" << G4endl;
    for(; itrMap != map->GetMap()->end(); itrMap++) {
      idx = itrMap->first;
      GetXYZ(idx, q);
      G4cout << "  (" << std::setw(3) << q[0] << ","
	     << std::setw(3) << q[1] << ","
	     << std::setw(3) << q[2] << ")     "
	     << *(itrMap->second) << G4endl;
      tot += *(itrMap->second);
    }
    G4cout << "       total :    " << tot << G4endl;
    */
  }
  //G4cout << G4endl;


}

void G4ScoringBox::Draw() {

  G4int nps = fMFD->GetNumberOfPrimitives(); 
  G4PSEnergyDeposit3D * ed3d;
  for(int i = 0; i < nps; i++) {
    G4VPrimitiveScorer * ps = fMFD->GetPrimitive(i);
    ed3d = dynamic_cast<G4PSEnergyDeposit3D*>(ps);
    if(ed3d != NULL) break;
  }

  G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager && ed3d) {
    
    // edep vectors
    std::vector<std::vector<std::vector<double> > > edep; // edep[X][Y][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[2]; z++) ez.push_back(0.);
    std::vector<std::vector<double> > eyz;
    for(int y = 0; y < fNSegment[1]; y++) eyz.push_back(ez);
    for(int x = 0; x < fNSegment[0]; x++) edep.push_back(eyz);

    std::vector<std::vector<double> > xyedep; // xyedep[X][Y]
    std::vector<double> ey;
    for(int y = 0; y < fNSegment[1]; y++) ey.push_back(0.);
    for(int x = 0; x < fNSegment[0]; x++) xyedep.push_back(ey);

    std::vector<std::vector<double> > yzedep; // yzedep[Y][Z]
    for(int y = 0; y < fNSegment[1]; y++) yzedep.push_back(ez);

    std::vector<std::vector<double> > xzedep; // xzedep[X][Z]
    for(int x = 0; x < fNSegment[0]; x++) xzedep.push_back(ez);


    std::map<G4int, G4double*> * map = fMap.find(ed3d->GetName())->second->GetMap();
    G4double xymax = 0., yzmax = 0., xzmax = 0.;
    G4int q[3];
    std::map<G4int, G4double*>::iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      GetXYZ(itr->first, q);

      xyedep[q[0]][q[1]] += *(itr->second);
      if(xymax < xyedep[q[0]][q[1]]) xymax = xyedep[q[0]][q[1]];

      yzedep[q[1]][q[2]] += *(itr->second);
      if(yzmax < yzedep[q[1]][q[2]]) yzmax = yzedep[q[1]][q[2]];

      xzedep[q[0]][q[2]] += *(itr->second);
      if(xzmax < xzedep[q[0]][q[1]]) xzmax = xzedep[q[1]][q[2]];
    }  
    
    G4Box box("hitsbox", fSize[0]/fNSegment[0]*0.99,
	               fSize[1]/fNSegment[1]*0.99,
                       fSize[2]/fNSegment[2]*0.99);
    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);

    // xy plane
    for(int x = 0; x < fNSegment[0]; x++) {
      for(int y = 0; y < fNSegment[1]; y++) {
        G4ThreeVector pos(GetReplicaPosition(x, y, 0) - fCenterPosition);
        G4ThreeVector pos2(GetReplicaPosition(x, y, fNSegment[2]-1) - fCenterPosition);
        G4Transform3D trans, trans2;
        if(fRotationMatrix) {
          trans = G4Transform3D(*fRotationMatrix, pos);
          trans2 = G4Transform3D(*fRotationMatrix, pos2);
        } else {
          trans = G4Translate3D(pos);
          trans2 = G4Translate3D(pos2);
        }
	G4double c[4];
	GetMapColor(xyedep[x][y]/xymax, c);
	att.SetColour(c[0], c[1], c[2]);//, c[3]);
	pVisManager->Draw(box, att, trans);
	pVisManager->Draw(box, att, trans2);

      }
    }

    // yz plane
    for(int y = 0; y < fNSegment[1]; y++) {
      for(int z = 0; z < fNSegment[2]; z++) {
        G4ThreeVector pos(GetReplicaPosition(0, y, z) - fCenterPosition);
        G4ThreeVector pos2(GetReplicaPosition(fNSegment[0]-1, y, z) - fCenterPosition);
        G4Transform3D trans, trans2;
        if(fRotationMatrix) {
          trans = G4Transform3D(*fRotationMatrix, pos);
          trans2 = G4Transform3D(*fRotationMatrix, pos2);
        } else {
          trans = G4Translate3D(pos);
          trans2 = G4Translate3D(pos2);
        }
	G4double c[4];
	GetMapColor(yzedep[y][z]/yzmax, c);
	att.SetColour(c[0], c[1], c[2]);//, c[3]);
	pVisManager->Draw(box, att, trans);
	pVisManager->Draw(box, att, trans2);

      }
    }

    // xz plane
    for(int x = 0; x < fNSegment[0]; x++) {
      for(int z = 0; z < fNSegment[2]; z++) {
        G4ThreeVector pos(GetReplicaPosition(x, 0, z) - fCenterPosition );
        G4ThreeVector pos2(GetReplicaPosition(x, fNSegment[1]-1, z) - fCenterPosition);
        G4Transform3D trans, trans2;
        if(fRotationMatrix) {
          trans = G4Transform3D(*fRotationMatrix, pos);
          trans2 = G4Transform3D(*fRotationMatrix, pos2);
        } else {
          trans = G4Translate3D(pos);
          trans2 = G4Translate3D(pos2);
        }
	G4double c[4];
	GetMapColor(xzedep[x][z]/xzmax, c);
	att.SetColour(c[0], c[1], c[2]);//, c[3]);
	pVisManager->Draw(box, att, trans);
	pVisManager->Draw(box, att, trans2);

      }
    }

  }
}

G4ThreeVector G4ScoringBox::GetReplicaPosition(G4int x, G4int y, G4int z) {
  G4ThreeVector width(fSize[0]/fNSegment[0], fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]);
  G4ThreeVector pos(-fSize[0] + 2*(x+0.5)*width.x(),
		    -fSize[1] + 2*(y+0.5)*width.y(),
		    -fSize[2] + 2*(z+0.5)*width.z());

  return pos;
}

/*
void G4ScoringBox::GetSegmentOrder(G4int segDir, G4int nseg[3], G4int segOrd[3], G4double segfact[3][3]) {

  if(segDir == 1) { // in x direction, it segments z -> y -> x directions by turns.
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

  } else if(segDir == -1 || segDir == 3) { // in z direction, it segments x -> y -> z directions by turns.
    segOrd[0] = 0;
    segOrd[1] = 1;
    segOrd[2] = 2;
    segfact[0][0] = G4double(nseg[0]);
    segfact[0][1] = 1.;
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
*/

void G4ScoringBox::GetXYZ(G4int index, G4int q[3]) const {

  q[0] = index/(fNSegment[2]*fNSegment[1]);
  q[1] = (index - q[0]*fNSegment[2]*fNSegment[1])/fNSegment[1];
  q[2] = index - q[1]*fNSegment[1] - q[0]*fNSegment[2]*fNSegment[1];
  /*
  if(fSegmentDirection == 3 || fSegmentDirection == -1) {
    q[2] = index/(fNSegment[0]*fNSegment[1]);
    q[1] = (index - q[2]*fNSegment[0]*fNSegment[1])/fNSegment[1];
    q[0] = index - q[1]*fNSegment[0] - q[2]*fNSegment[0]*fNSegment[1];

  } else if(fSegmentDirection == 2) {
    q[1] = index/(fNSegment[0]*fNSegment[2]);
    q[2] = (index - q[1]*fNSegment[0]*fNSegment[2])/fNSegment[2];
    q[0] = index - q[2]*fNSegment[0] - q[1]*fNSegment[0]*fNSegment[2];
    
  } else if(fSegmentDirection == 1) {
    q[0] = index/(fNSegment[2]*fNSegment[1]);
    q[1] = (index - q[0]*fNSegment[2]*fNSegment[1])/fNSegment[1];
    q[2] = index - q[1]*fNSegment[2] - q[0]*fNSegment[2]*fNSegment[1];
    
  }
  */

  //G4cout << "GetXYZ: " << index << ": "
  //<< q[0] << ", " << q[1] << ", " << q[2] << G4endl;
}

void G4ScoringBox::GetMapColor(G4double value, G4double color[4]) {
  color[0] = color[1] = color[2] = 0.;
  color[4] = 1.;

  if(value > 0.8) {
    color[0] = 1.;
    color[1] = 1.-(value-0.7)/0.3;
  } else if(value > 0.6) {
    color[0] = (value-0.6)/0.2;
    color[1] = 1.;
  } else if(value > 0.4) {
    color[1] = 1.;
    color[2] = 1.-(value-0.4)/0.2;
  } else if(value > 0.2) {
    color[2] = 1.;
    color[0] = color[1] = 1. - (value-0.2)/0.2;
  } else {
    color[2] = 1.;
    color[0] = color[1] = 1. - value/0.2;
  }
}
