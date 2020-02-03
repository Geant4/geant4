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
//

#include "G4ScoringCylinder.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4VScoreColorMap.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"
#include "G4PSNofStep.hh"
#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4ScoringCylinder::G4ScoringCylinder(G4String wName)
  :G4VScoringMesh(wName)
{
  fShape = MeshShape::cylinder;

  fDivisionAxisNames[0] = "Z";
  fDivisionAxisNames[1] = "PHI";
  fDivisionAxisNames[2] = "R";
}

G4ScoringCylinder::~G4ScoringCylinder()
{;}

void G4ScoringCylinder::SetupGeometry(G4VPhysicalVolume * fWorldPhys) {

  if(verboseLevel > 9) G4cout << "G4ScoringCylinder::SetupGeometry() ..." << G4endl;

  // World
  G4VPhysicalVolume * scoringWorld = fWorldPhys;
  G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  if(verboseLevel > 9) G4cout << fWorldName << G4endl;
  G4String tubsName = fWorldName+"_mesh";

  if(verboseLevel > 9) G4cout << "R max., Dz =: " << fSize[0] << ", " << fSize[1] << G4endl;
  G4VSolid * tubsSolid = new G4Tubs(tubsName+"0", // name
				    0.,           // R min
				    fSize[0],     // R max
				    fSize[1],     // Dz
				    0.,           // starting phi
                                    twopi*rad);   // segment phi
  G4LogicalVolume *  tubsLogical = new G4LogicalVolume(tubsSolid, 0, tubsName);
  new G4PVPlacement(fRotationMatrix, fCenterPosition,
		    tubsLogical, tubsName+"0", worldLogical, false, 0);

  if(verboseLevel > 9) G4cout << " # of segments : r, phi, z =: "
	 << fNSegment[IR] << ", " << fNSegment[IPHI] << ", " << fNSegment[IZ] << G4endl;

  G4String layerName[2] = {tubsName + "1",  tubsName + "2"};
  G4VSolid * layerSolid[2]; 
  G4LogicalVolume * layerLogical[2];

  //-- fisrt nested layer (replicated along z direction)
  if(verboseLevel > 9) G4cout << "layer 1 :" << G4endl;
  layerSolid[0] = new G4Tubs(layerName[0],           // name
			     0.,                     // inner radius
			     fSize[0],               // outer radius
			     fSize[1]/fNSegment[IZ], // half len. in z
			     0.,                     // starting phi angle
			     twopi*rad);             // delta angle of the segment
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], 0, layerName[0]);
  if(fNSegment[IZ] > 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replicate along z direction" << G4endl;
    if(G4ScoringManager::GetReplicaLevel()>0) {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replica" << G4endl;
      new G4PVReplica(layerName[0], layerLogical[0], tubsLogical, kZAxis, fNSegment[IZ], 2.*fSize[1]/fNSegment[IZ]);
    } else {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Division" << G4endl;
      new G4PVDivision(layerName[0], layerLogical[0], tubsLogical, kZAxis, fNSegment[IZ], 0.);
    }
  } else if(fNSegment[IZ] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[0], layerName[0], tubsLogical, false, 0);
  } else {
    G4cerr << "G4ScoringCylinder::SetupGeometry() : invalid parameter ("
	   << fNSegment[IZ] << ") "
	   << "in placement of the first nested layer." << G4endl;
  }

  // second nested layer (replicated along phi direction)
  if(verboseLevel > 9) G4cout << "layer 2 :" << G4endl;
  layerSolid[1] = new G4Tubs(layerName[1],
			     0.,
			     fSize[0],
			     fSize[1]/fNSegment[IZ],
			     0.,
                             twopi*rad/fNSegment[IPHI]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], 0, layerName[1]);
  if(fNSegment[IPHI] > 1)  {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replicate along phi direction" << G4endl;
    if(G4ScoringManager::GetReplicaLevel()>1) {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replica" << G4endl;
      new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], kPhi,
		      fNSegment[IPHI], twopi*rad/fNSegment[IPHI]);
    } else {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Division" << G4endl;
      new G4PVDivision(layerName[1], layerLogical[1], layerLogical[0], kPhi, fNSegment[IPHI], 0.);
    }
  } else if(fNSegment[IPHI] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), layerLogical[1], layerName[1], layerLogical[0], false, 0);
  } else
    G4cerr << "ERROR : G4ScoringCylinder::SetupGeometry() : invalid parameter ("
	   << fNSegment[IPHI] << ") "
	   << "in placement of the second nested layer." << G4endl;

  // mesh elements
  if(verboseLevel > 9) G4cout << "mesh elements :" << G4endl;
  G4String elementName = tubsName +"3";
  G4VSolid * elementSolid = new G4Tubs(elementName,
				       0.,
				       fSize[0]/fNSegment[IR],
				       fSize[1]/fNSegment[IZ],
				       0.,
                                       twopi*rad/fNSegment[IPHI]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, 0, elementName);
  if(fNSegment[IR] > 1) {

    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replicate along r direction" << G4endl;

    if(G4ScoringManager::GetReplicaLevel()>2) {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Replica" << G4endl;
      new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], kRho,
			fNSegment[IR], fSize[0]/fNSegment[IR]);
    } else {
      if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Division" << G4endl;
      new G4PVDivision(elementName, fMeshElementLogical, layerLogical[1], kRho, fNSegment[IR], 0.);
    }
  } else if(fNSegment[IR] == 1) {
    if(verboseLevel > 9) G4cout << "G4ScoringCylinder::Construct() : Placement" << G4endl;
    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), fMeshElementLogical, elementName, layerLogical[1], false, 0);
  } else {
    G4cerr << "G4ScoringCylinder::SetupGeometry() : "
	   << "invalid parameter (" << fNSegment[IR] << ") "
	   << "in mesh element placement." << G4endl;
  }

  // set the sensitive detector
  fMeshElementLogical->SetSensitiveDetector(fMFD);
  

  // vis. attributes
  G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5));
  visatt->SetVisibility(true);
  layerLogical[0]->SetVisAttributes(visatt);
  layerLogical[1]->SetVisAttributes(visatt);
  visatt = new G4VisAttributes(G4Colour(.5,.5,.5,0.01));
  //visatt->SetForceSolid(true);
  fMeshElementLogical->SetVisAttributes(visatt);
}

void G4ScoringCylinder::List() const {
  G4cout << "G4ScoringCylinder : " << fWorldName << " --- Shape: Cylindrical mesh" << G4endl;

  G4cout << " Size (R, Dz): ("
	 << fSize[0]/cm << ", "
	 << fSize[1]/cm << ") [cm]"
	 << G4endl;

  G4VScoringMesh::List();
}


void G4ScoringCylinder::Draw(RunScore * map,
                             G4VScoreColorMap* colorMap, G4int axflg) {

  G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager) {
    
    // cell vectors
    std::vector<double> ephi;
    for(int phi = 0; phi < fNSegment[IPHI]; phi++) ephi.push_back(0.);
    //-
    std::vector<std::vector<double> > zphicell; // zphicell[Z][PHI]
    for(int z = 0; z < fNSegment[IZ]; z++) zphicell.push_back(ephi);
    //-
    std::vector<std::vector<double> > rphicell; // rphicell[R][PHI]
    for(int r = 0; r < fNSegment[IR]; r++) rphicell.push_back(ephi);
    
    // projections
    G4int q[3];
    std::map<G4int, G4StatDouble*>::iterator itr = map->GetMap()->begin();
    for(; itr != map->GetMap()->end(); itr++) {
      if(itr->first < 0) {
        G4cout << itr->first << G4endl;
        continue;
      }
      GetRZPhi(itr->first, q);
      
      zphicell[q[IZ]][q[IPHI]] += (itr->second->sum_wx())/fDrawUnitValue;
      rphicell[q[IR]][q[IPHI]] += (itr->second->sum_wx())/fDrawUnitValue;
    }
    
    // search min./max. values
    G4double zphimin = DBL_MAX, rphimin = DBL_MAX;
    G4double zphimax = 0., rphimax = 0.;
    for(int iphi = 0; iphi < fNSegment[IPHI]; iphi++) {
      for(int iz = 0; iz < fNSegment[IZ]; iz++) {
        if(zphimin > zphicell[iz][iphi]) zphimin = zphicell[iz][iphi];
        if(zphimax < zphicell[iz][iphi]) zphimax = zphicell[iz][iphi];
      }
      for(int ir = 0; ir < fNSegment[IR]; ir++) {
        if(rphimin > rphicell[ir][iphi]) rphimin = rphicell[ir][iphi];
        if(rphimax < rphicell[ir][iphi]) rphimax = rphicell[ir][iphi];
      }
    }
    
    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);
    
    
    G4Scale3D scale;
    if(axflg/100==1) {
      // rz plane
    }
    axflg = axflg%100;
    if(axflg/10==1) {
      pVisManager->BeginDraw();
      
      // z-phi plane
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(zphimin, zphimax); }
      
      G4double zhalf = fSize[1]/fNSegment[IZ];
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        for(int z = 0; z < fNSegment[IZ]; z++) {
          //-
          G4double angle = twopi/fNSegment[IPHI]*phi;
          G4double dphi = twopi/fNSegment[IPHI];
          G4Tubs cylinder("z-phi",                    // name
                          fSize[0]*0.99, fSize[0],  // inner radius, outer radius
                          zhalf,                      // half length in z
                          angle, dphi*0.99999);       // starting phi angle, delta angle
          //-
          G4ThreeVector zpos(0., 0., -fSize[1] + fSize[1]/fNSegment[IZ]*(1 + 2.*z));
          G4Transform3D trans;
          if(fRotationMatrix) {
            trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zpos);
            trans = G4Translate3D(fCenterPosition)*trans;
          } else {
            trans = G4Translate3D(zpos)*G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(zphicell[z][phi], c);
          att.SetColour(c[0], c[1], c[2]);//, c[3]);
          //-
          G4Polyhedron * poly = cylinder.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
      pVisManager->EndDraw();
    }
    axflg = axflg%10;
    if(axflg==1) {
      pVisManager->BeginDraw();
      
      // r-phi plane
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(rphimin, rphimax); }
      
      G4double rsize = fSize[0]/fNSegment[IR];
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        for(int r = 0; r < fNSegment[IR]; r++) {
          
          G4double rs[2] = {rsize*r, rsize*(r+1)};
          G4double angle = twopi/fNSegment[IPHI]*phi;
          G4double dphi = twopi/fNSegment[IPHI];
          G4Tubs cylindern("z-phi", rs[0], rs[1], 0.001,
                          angle, dphi*0.99999);
          G4Tubs cylinderp = cylindern;
          
          G4ThreeVector zposn(0., 0., -fSize[1]);
          G4ThreeVector zposp(0., 0.,  fSize[1]);
          G4Transform3D transn, transp;
          if(fRotationMatrix) {
            transn = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zposn);
            transn = G4Translate3D(fCenterPosition)*transn;
            transp = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zposp);
            transp = G4Translate3D(fCenterPosition)*transp;
          } else {
            transn = G4Translate3D(zposn)*G4Translate3D(fCenterPosition);
            transp = G4Translate3D(zposp)*G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(rphicell[r][phi], c);
          att.SetColour(c[0], c[1], c[2]);//, c[3]);

          G4Polyhedron * polyn = cylindern.GetPolyhedron();
          polyn->Transform(transn);
          polyn->SetVisAttributes(att);
          pVisManager->Draw(*polyn);

          G4Polyhedron * polyp = cylinderp.GetPolyhedron();
          polyp->Transform(transp);
          polyp->SetVisAttributes(att);
          pVisManager->Draw(*polyp);
        }
      }
      
      pVisManager->EndDraw();
    }
    
    colorMap->SetPSUnit(fDrawUnit);
    colorMap->SetPSName(fDrawPSName);
    colorMap->DrawColorChart();
    
  }
}

void G4ScoringCylinder::DrawColumn(RunScore * map, G4VScoreColorMap* colorMap, 
				   G4int idxProj, G4int idxColumn) 
{
  G4int projAxis = 0;
  switch(idxProj) {
    case 0:
      projAxis = IR;
      break;
    case 1:
      projAxis = IZ;
      break;
    case 2:
      projAxis = IPHI;
      break;
  }
  
  if(idxColumn<0 || idxColumn>=fNSegment[projAxis])
  {
    G4cerr << "Warning : Column number " << idxColumn << " is out of scoring mesh [0," << fNSegment[projAxis]-1 <<
    "]. Method ignored." << G4endl;
    return;
  }
  G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager) {
    
    // cell vectors
    std::vector<std::vector<std::vector<double> > > cell; // cell[R][Z][PHI]
    std::vector<double> ephi;
    for(int phi = 0; phi < fNSegment[IPHI]; phi++) ephi.push_back(0.);
    std::vector<std::vector<double> > ezphi;
    for(int z = 0; z < fNSegment[IZ]; z++) ezphi.push_back(ephi);
    for(int r = 0; r < fNSegment[IR]; r++) cell.push_back(ezphi);
    
    std::vector<std::vector<double> > rzcell; // rzcell[R][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[IZ]; z++) ez.push_back(0.);
    for(int r = 0; r < fNSegment[IR]; r++) rzcell.push_back(ez);
    
    std::vector<std::vector<double> > zphicell; // zphicell[Z][PHI]
    for(int z = 0; z < fNSegment[IZ]; z++) zphicell.push_back(ephi);
    
    std::vector<std::vector<double> > rphicell; // rphicell[R][PHI]
    for(int r = 0; r < fNSegment[IR]; r++) rphicell.push_back(ephi);
    
    // projections
    G4int q[3];
    std::map<G4int,G4StatDouble*>::iterator itr = map->GetMap()->begin();
    for(; itr != map->GetMap()->end(); itr++) {
      if(itr->first < 0) {
        G4cout << itr->first << G4endl;
        continue;
      }
      GetRZPhi(itr->first, q);
      
      if(projAxis == IR && q[IR] == idxColumn) { // zphi plane
        zphicell[q[IZ]][q[IPHI]] += (itr->second->sum_wx())/fDrawUnitValue;
      }
      if(projAxis == IZ && q[IZ] == idxColumn) { // rphi plane
        rphicell[q[IR]][q[IPHI]] += (itr->second->sum_wx())/fDrawUnitValue;
      }
      if(projAxis == IPHI && q[IPHI] == idxColumn) { // rz plane
        rzcell[q[IR]][q[IZ]] += (itr->second->sum_wx())/fDrawUnitValue;
      }
    }
    
    // search min./max. values
    G4double rzmin = DBL_MAX, zphimin = DBL_MAX, rphimin = DBL_MAX;
    G4double rzmax = 0., zphimax = 0., rphimax = 0.;
    for(int r = 0; r < fNSegment[IR]; r++) {
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        if(rphimin > rphicell[r][phi]) rphimin = rphicell[r][phi];
        if(rphimax < rphicell[r][phi]) rphimax = rphicell[r][phi];
      }
      for(int z = 0; z < fNSegment[IZ]; z++) {
        if(rzmin > rzcell[r][z]) rzmin = rzcell[r][z];
        if(rzmax < rzcell[r][z]) rzmax = rzcell[r][z];
      }
    }
    for(int z = 0; z < fNSegment[IZ]; z++) {
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        if(zphimin > zphicell[z][phi]) zphimin = zphicell[z][phi];
        if(zphimax < zphicell[z][phi]) zphimax = zphicell[z][phi];
      }
    }
    
    
    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);
    
    pVisManager->BeginDraw();
    
    G4Scale3D scale;
    // z-phi plane
    if(projAxis == IR) {
      
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(zphimin,zphimax); }
      
      G4double zhalf = fSize[1]/fNSegment[IZ];
      G4double rsize[2] = {fSize[0]/fNSegment[IR]*idxColumn,
        fSize[0]/fNSegment[IR]*(idxColumn+1)};
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        for(int z = 0; z < fNSegment[IZ]; z++) {
          
          G4double angle = twopi/fNSegment[IPHI]*phi*radian;
          G4double dphi = twopi/fNSegment[IPHI]*radian;
          G4Tubs cylinder("z-phi", rsize[0], rsize[1], zhalf,
                          angle, dphi*0.99999);
          
          G4ThreeVector zpos(0., 0., -fSize[1] + fSize[1]/fNSegment[IZ]*(1 + 2.*z));
          G4Transform3D trans;
          if(fRotationMatrix) {
            trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zpos);
            trans = G4Translate3D(fCenterPosition)*trans;
          } else {
            trans = G4Translate3D(zpos)*G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(zphicell[z][phi], c);
          att.SetColour(c[0], c[1], c[2]);//, c[3]);
          
          G4Polyhedron * poly = cylinder.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
      
      // r-phi plane
    } else if(projAxis == IZ) {
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(rphimin,rphimax); }
      
      G4double rsize = fSize[0]/fNSegment[IR];
      for(int phi = 0; phi < fNSegment[IPHI]; phi++) {
        for(int r = 0; r < fNSegment[IR]; r++) {
          
          G4double rs[2] = {rsize*r, rsize*(r+1)};
          G4double angle = twopi/fNSegment[IPHI]*phi*radian;
          G4double dz = fSize[1]/fNSegment[IZ];
          G4double dphi = twopi/fNSegment[IPHI]*radian;
          G4Tubs cylinder("r-phi", rs[0], rs[1], dz,
                          angle, dphi*0.99999);
          G4ThreeVector zpos(0., 0.,
                             -fSize[1]+fSize[1]/fNSegment[IZ]*(idxColumn*2+1));
          G4Transform3D trans;
          if(fRotationMatrix) {
            trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zpos);
            trans = G4Translate3D(fCenterPosition)*trans;
          } else {
            trans = G4Translate3D(zpos)*G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(rphicell[r][phi], c);
          att.SetColour(c[0], c[1], c[2]);//, c[3]);

          G4Polyhedron * poly = cylinder.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
      
      // r-z plane
    } else if(projAxis == IPHI) {
      if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(rzmin,rzmax); }
      
      G4double rsize = fSize[0]/fNSegment[IR];
      G4double zhalf = fSize[1]/fNSegment[IZ];
      G4double angle = twopi/fNSegment[IPHI]*idxColumn*radian;
      G4double dphi = twopi/fNSegment[IPHI]*radian;
      for(int z = 0; z < fNSegment[IZ]; z++) {
        for(int r = 0; r < fNSegment[IR]; r++) {
          
          G4double rs[2] = {rsize*r, rsize*(r+1)};
          G4Tubs cylinder("z-phi", rs[0], rs[1], zhalf,
                          angle, dphi);
          
          G4ThreeVector zpos(0., 0.,
                             -fSize[1]+fSize[1]/fNSegment[IZ]*(2.*z+1));
          G4Transform3D trans;
          if(fRotationMatrix) {
            trans = G4Rotate3D(*fRotationMatrix).inverse()*G4Translate3D(zpos);
            trans = G4Translate3D(fCenterPosition)*trans;
          } else {
            trans = G4Translate3D(zpos)*G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(rzcell[r][z], c);
          att.SetColour(c[0], c[1], c[2]);//, c[3]);

          G4Polyhedron * poly = cylinder.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
    }
    pVisManager->EndDraw();
  }
  
  colorMap->SetPSUnit(fDrawUnit);
  colorMap->SetPSName(fDrawPSName);
  colorMap->DrawColorChart();
  
}

void G4ScoringCylinder::GetRZPhi(G4int index, G4int q[3]) const {
  // index = k + j * k-size + i * jk-plane-size

  // nested : z -> phi -> r
  G4int i = IZ;
  G4int j = IPHI;
  G4int k = IR;
  G4int jk = fNSegment[j]*fNSegment[k];
  q[i] = index/jk;
  q[j] = (index - q[i]*jk)/fNSegment[k];
  q[k] = index - q[j]*fNSegment[k] - q[i]*jk;
}
