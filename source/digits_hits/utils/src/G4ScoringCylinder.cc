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
// $Id: G4ScoringCylinder.cc,v 1.2 2008-03-09 12:28:40 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScoringCylinder.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4PVParameterised.hh"
#include "G4VisAttributes.hh"
//#include "G4ScoringCylinderParameterisation.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4ScoringManager.hh"


G4ScoringCylinder::G4ScoringCylinder(G4String wName)
  :G4VScoringMesh(wName), fSegmentDirection(-1), 
   fMeshElementLogical(0)
{
  fShape = cylinderMesh;
}

G4ScoringCylinder::~G4ScoringCylinder()
{
}

void G4ScoringCylinder::Construct(G4VPhysicalVolume* fWorldPhys)
{
  if(fConstructed) {

    if(verboseLevel > 0) 
      G4cout << fWorldPhys->GetName() << " --- All quantities are reset." << G4endl;
    ResetScore();

  } else {
    fConstructed = true;
    SetupGeometry(fWorldPhys);
  }
}



void G4ScoringCylinder::SetupGeometry(G4VPhysicalVolume * fWorldPhys) {

  if(verboseLevel > 9) G4cout << "G4ScoringCylinder::SetupGeometry() ..." << G4endl;

  // World
  G4VPhysicalVolume * scoringWorld = fWorldPhys;
  G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  if(verboseLevel > 9) G4cout << fWorldName << G4endl;
  G4String tubsName = fWorldName;

  if(verboseLevel > 9) G4cout << fSize[0] << ", " << fSize[1] << G4endl;
  G4VSolid * tubsSolid = new G4Tubs(tubsName+"0", 0., fSize[0], fSize[1],
				    0., twopi*rad);
  G4LogicalVolume *  tubsLogical = new G4LogicalVolume(tubsSolid, 0, tubsName);
  new G4PVPlacement(fRotationMatrix, fCenterPosition,
		    tubsLogical, tubsName+"0", worldLogical, false, 0);


  G4String layerName[2] = {tubsName + "1",  tubsName + "2"};
  G4VSolid * layerSolid[2]; 
  G4LogicalVolume * layerLogical[2];

  // fisrt nested layer (replicated along r direction)
  if(verboseLevel > 9) G4cout << "layer 1 :" << G4endl;
  layerSolid[0] = new G4Tubs(layerName[0],
			     0.,
			     fSize[0]/fNSegment[0],
			     fSize[1],
			     0., twopi);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  if(fNSegment[0] > 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replicate along r direction" << G4endl;
    G4double r = fSize[0]/fNSegment[0];
    if(G4ScoringManager::GetReplicaLevel()>0) {

      new G4PVReplica(layerName[0], layerLogical[0], tubsLogical, kRho,
		      fNSegment[0], r, 0.);
    } else {
      new G4PVDivision(layerName[0], layerLogical[0], tubsLogical, kRho,
		       fNSegment[0], r);
    }
  } else if(fNSegment[0] == 1) {
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[0], layerName[0], tubsLogical, false, 0);
  } else {
    G4cerr << "G4ScoringCylinder::SetupGeometry() : invalid parameter ("
	   << fNSegment[0] << ") "
	   << "in placement of the first nested layer." << G4endl;
  }

  if(verboseLevel > 9) {
    G4cout << fSize[0] << ", "
	   << fSize[1] 
	   << G4endl;
    G4cout << layerName[0] << ": kRho, "
	   << fNSegment[0] << ", "
	   << fSize[0]/fNSegment[0] << G4endl;
  }

  // second nested layer (replicated along z direction)
  if(verboseLevel > 9) G4cout << "layer 2 :" << G4endl;
  layerSolid[1] = new G4Tubs(layerName[1],
			     0.,
			     fSize[0],///fNSegment[0],
			     fSize[1]/fNSegment[1],
			     0., twopi);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  if(fNSegment[1] > 1)  {
    G4double width = fSize[1]/fNSegment[1]*2.;
    if(G4ScoringManager::GetReplicaLevel()>1) {

      new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], kZAxis,
		      fNSegment[1], width);
    } else {
      new G4PVDivision(layerName[1], layerLogical[1], layerLogical[0], kZAxis,
		       fNSegment[1], width);
    }
  } else if(fNSegment[1] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[1], layerName[1], layerLogical[0], false, 0);
  } else
    G4cerr << "ERROR : G4ScoringCylinder::SetupGeometry() : invalid parameter ("
	   << fNSegment[1] << ") "
	   << "in placement of the second nested layer." << G4endl;

  if(verboseLevel > 9) {
    G4cout << fSize[0]/fNSegment[0] << ", "
	   << fSize[1]/fNSegment[1] << G4endl;
    G4cout << layerName[1] << ": kZAxis, "
	   << fNSegment[1] << ", "
	   << fSize[1]/fNSegment[1] << G4endl;
  }


  // mesh elements
  if(verboseLevel > 9) G4cout << "mesh elements :" << G4endl;
  G4String elementName = tubsName +"3";
  G4VSolid * elementSolid = new G4Tubs(elementName,
				       0.,
				       fSize[0],//fNSegment[0],
				       fSize[1]/fNSegment[1],
				       0., twopi/fNSegment[2]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, 0, elementName);
  if(fNSegment[2] > 1) {

    /*
    if(fSegmentPositions.size() > 0) {
      G4double motherDims[3] ={fSize[0]/fsegParam[2][0],
			       fSize[1]/fsegParam[2][1],
			       fSize[2]/fsegParam[2][2]};
      G4int nelement = fSegmentPositions.size() + 1;
      //G4ScoringCylinderParameterisation * param =
      G4VPVParameterisation * param =
	new G4ScoringCylinderParameterisation(axis[2], motherDims, fSegmentPositions);
      new G4PVParameterised(elementName,
			    fMeshElementLogical,
			    layerLogical[1],
			    axis[2],
			    nelement,
			    param);
    } else {
    */
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replicate along phi direction" << G4endl;

    G4double angle = twopi*rad/fNSegment[2];
    //if(G4ScoringManager::GetReplicaLevel()>2) {
    if(1) {
      new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], kPhi,
			fNSegment[2], angle, 0.);
    } else {
      new G4PVDivision(elementName, fMeshElementLogical, layerLogical[1], kPhi,
			 fNSegment[2], angle);
    }
    //}
  } else if(fNSegment[2] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), fMeshElementLogical, elementName, layerLogical[1], false, 0);
  } else {
    G4cerr << "G4ScoringCylinder::SetupGeometry() : "
	   << "invalid parameter (" << fNSegment[2] << ") "
	   << "in mesh element placement." << G4endl;
  }

  if(verboseLevel > 9) {
    G4cout << fSize[0]/fNSegment[0] << ", "
	   << fSize[1]/fNSegment[1] << G4endl;
    G4cout << elementName << ": kPhi, "
	   << fNSegment[2] << G4endl;
  }

  // set the sensitive detector
  fMeshElementLogical->SetSensitiveDetector(fMFD);
  

  // vis. attributes
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5,0.1));
  visatt->SetVisibility(true);
  layerLogical[0]->SetVisAttributes(visatt);
  layerLogical[1]->SetVisAttributes(visatt);
  visatt = new G4VisAttributes(G4Colour(.5,.5,.5,0.01));
  visatt->SetForceSolid(true);
  fMeshElementLogical->SetVisAttributes(visatt);
}

void G4ScoringCylinder::List() const {
  G4cout << "G4ScoringCylinder : " << fWorldName << " --- Shape: Cylindrical mesh" << G4endl;

  G4VScoringMesh::List();
}


void G4ScoringCylinder::Draw(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap, G4int axflg) {

  /*
  G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager) {
    
    // cell vectors
    std::vector<std::vector<std::vector<double> > > cell; // cell[X][Y][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[2]; z++) ez.push_back(0.);
    std::vector<std::vector<double> > eyz;
    for(int y = 0; y < fNSegment[1]; y++) eyz.push_back(ez);
    for(int x = 0; x < fNSegment[0]; x++) cell.push_back(eyz);

    std::vector<std::vector<double> > xycell; // xycell[X][Y]
    std::vector<double> ey;
    for(int y = 0; y < fNSegment[1]; y++) ey.push_back(0.);
    for(int x = 0; x < fNSegment[0]; x++) xycell.push_back(ey);

    std::vector<std::vector<double> > yzcell; // yzcell[Y][Z]
    for(int y = 0; y < fNSegment[1]; y++) yzcell.push_back(ez);

    std::vector<std::vector<double> > xzcell; // xzcell[X][Z]
    for(int x = 0; x < fNSegment[0]; x++) xzcell.push_back(ez);

    G4double xymax = 0., yzmax = 0., xzmax = 0.;
    G4int q[3];
    std::map<G4int, G4double*>::iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      GetXYZ(itr->first, q);

      xycell[q[0]][q[1]] += *(itr->second);
      if(xymax < xycell[q[0]][q[1]]) xymax = xycell[q[0]][q[1]];

      yzcell[q[1]][q[2]] += *(itr->second);
      if(yzmax < yzcell[q[1]][q[2]]) yzmax = yzcell[q[1]][q[2]];

      xzcell[q[0]][q[2]] += *(itr->second);
      if(xzmax < xzcell[q[0]][q[2]]) xzmax = xzcell[q[0]][q[2]];
    }  
    
    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);


    G4Scale3D scale;
    if(axflg/100==1) {
      // xy plane
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,xymax); }
      G4ThreeVector zhalf(0., 0., fSize[2]/fNSegment[2]*0.98);
      G4Box xyplate("xy", fSize[0]/fNSegment[0], fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]*0.01);
      for(int x = 0; x < fNSegment[0]; x++) {
	for(int y = 0; y < fNSegment[1]; y++) {
	  G4ThreeVector pos(GetReplicaPosition(x, y, 0) - zhalf);
	  G4ThreeVector pos2(GetReplicaPosition(x, y, fNSegment[2]-1) + zhalf);
	  G4Transform3D trans, trans2;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	    trans2 = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos2);
	    trans2 = G4Translate3D(fCenterPosition)*trans2;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	    trans2 = G4Translate3D(pos2)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(xycell[x][y], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(xyplate, att, trans);
	  pVisManager->Draw(xyplate, att, trans2);

	}
      }
    }
    axflg = axflg%100;
    if(axflg/10==1) {
      // yz plane
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,yzmax); }
      G4ThreeVector xhalf(fSize[0]/fNSegment[0]*0.98, 0., 0.);
      G4Box yzplate("yz", fSize[0]/fNSegment[0]*0.01, fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]);
      for(int y = 0; y < fNSegment[1]; y++) {
	for(int z = 0; z < fNSegment[2]; z++) {
	  G4ThreeVector pos(GetReplicaPosition(0, y, z) - xhalf);
	  G4ThreeVector pos2(GetReplicaPosition(fNSegment[0]-1, y, z) + xhalf);
	  G4Transform3D trans, trans2;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	    trans2 = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos2);
	    trans2 = G4Translate3D(fCenterPosition)*trans2;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	    trans2 = G4Translate3D(pos2)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(yzcell[y][z], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(yzplate, att, trans);
	  pVisManager->Draw(yzplate, att, trans2);

	}
      }
    }
    axflg = axflg%10;
    if(axflg==1) {
      // xz plane
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,xzmax); }
      G4ThreeVector yhalf(0., fSize[1]/fNSegment[1]*0.98, 0.);
      G4Box xzplate("xz", fSize[0]/fNSegment[0], fSize[1]/fNSegment[1]*0.01, fSize[2]/fNSegment[2]);
      for(int x = 0; x < fNSegment[0]; x++) {
	for(int z = 0; z < fNSegment[2]; z++) {
	  G4ThreeVector pos(GetReplicaPosition(x, 0, z) - yhalf);
	  G4ThreeVector pos2(GetReplicaPosition(x, fNSegment[1]-1, z) + yhalf);
	  G4Transform3D trans, trans2;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	    trans2 = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos2);
	    trans2 = G4Translate3D(fCenterPosition)*trans2;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	    trans2 = G4Translate3D(pos2)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(xzcell[x][z], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(xzplate, att, trans);
	  pVisManager->Draw(xzplate, att, trans2);

	}
      }
    }
  }
  colorMap->DrawColorChart();
  */
}

void G4ScoringCylinder::DrawColumn(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap, 
				   G4int idxProj, G4int idxColumn) 
{

  /*
  if(idxColumn<0 || idxColumn>=fNSegment[idxProj])
  {
    G4cerr << "ERROR : Column number " << idxColumn << " is out of scoring mesh [0," << fNSegment[idxProj]-1 <<
    "]. Method ignored." << G4endl;
    return;
  }
  G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager) {
    
    // cell vectors
    std::vector<std::vector<std::vector<double> > > cell; // cell[X][Y][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[2]; z++) ez.push_back(0.);
    std::vector<std::vector<double> > eyz;
    for(int y = 0; y < fNSegment[1]; y++) eyz.push_back(ez);
    for(int x = 0; x < fNSegment[0]; x++) cell.push_back(eyz);

    std::vector<std::vector<double> > xycell; // xycell[X][Y]
    std::vector<double> ey;
    for(int y = 0; y < fNSegment[1]; y++) ey.push_back(0.);
    for(int x = 0; x < fNSegment[0]; x++) xycell.push_back(ey);

    std::vector<std::vector<double> > yzcell; // yzcell[Y][Z]
    for(int y = 0; y < fNSegment[1]; y++) yzcell.push_back(ez);

    std::vector<std::vector<double> > xzcell; // xzcell[X][Z]
    for(int x = 0; x < fNSegment[0]; x++) xzcell.push_back(ez);

    G4double xymax = 0., yzmax = 0., xzmax = 0.;
    G4int q[3];
    std::map<G4int, G4double*>::iterator itr = map->begin();
    for(; itr != map->end(); itr++) {
      GetXYZ(itr->first, q);

      if(idxProj == 0 && q[2] == idxColumn) { // xy plane
	xycell[q[0]][q[1]] += *(itr->second);
	if(xymax < xycell[q[0]][q[1]]) xymax = xycell[q[0]][q[1]];
      }
      if(idxProj == 1 && q[0] == idxColumn) { // yz plane
	yzcell[q[1]][q[2]] += *(itr->second);
	if(yzmax < yzcell[q[1]][q[2]]) yzmax = yzcell[q[1]][q[2]];
      }
      if(idxProj == 2 && q[1] == idxColumn) { // zx plane
	xzcell[q[0]][q[2]] += *(itr->second);
	if(xzmax < xzcell[q[0]][q[2]]) xzmax = xzcell[q[0]][q[2]];
      }
    }  
    
    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);


    G4Scale3D scale;
    // xy plane
    if(idxProj == 0) {
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,xymax); }
      G4Box xyplate("xy", fSize[0]/fNSegment[0], fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]);
      for(int x = 0; x < fNSegment[0]; x++) {
	for(int y = 0; y < fNSegment[1]; y++) {
	  G4ThreeVector pos(GetReplicaPosition(x, y, idxColumn));
	  G4Transform3D trans;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(xycell[x][y], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(xyplate, att, trans);

	}
      }
    } else
    // yz plane
    if(idxProj == 1) {
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,yzmax); }
      G4Box yzplate("yz", fSize[0]/fNSegment[0], fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]);
      for(int y = 0; y < fNSegment[1]; y++) {
	for(int z = 0; z < fNSegment[2]; z++) {
	  G4ThreeVector pos(GetReplicaPosition(idxColumn, y, z));
	  G4Transform3D trans;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(yzcell[y][z], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(yzplate, att, trans);
	}
      }
    } else
    // xz plane
      if(idxProj == 2) {
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(0.,xzmax);}
      G4Box xzplate("xz", fSize[0]/fNSegment[0], fSize[1]/fNSegment[1], fSize[2]/fNSegment[2]);
      for(int x = 0; x < fNSegment[0]; x++) {
	for(int z = 0; z < fNSegment[2]; z++) {
	  G4ThreeVector pos(GetReplicaPosition(x, idxColumn, z));
	  G4Transform3D trans;
	  if(fRotationMatrix) {
	    trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(pos);
	    trans = G4Translate3D(fCenterPosition)*trans;
	  } else {
	    trans = G4Translate3D(pos)*G4Translate3D(fCenterPosition);
	  }
	  G4double c[4];
	  colorMap->GetMapColor(xzcell[x][z], c);
	  att.SetColour(c[0], c[1], c[2]);//, c[3]);
	  pVisManager->Draw(xzplate, att, trans);
	}
      }
    }
  }

  colorMap->DrawColorChart();

  */
}


