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
// $Id: TiaraIronShieldB.cc,v 1.4 2006-06-29 15:45:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraIronShieldB.hh"
#include "TiaraDimensions.hh"
#include "TiaraMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"

/*
TiaraIronShieldB::
TiaraIronShieldB(TiaraMaterials &mfac,
		 const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetExperiment - 
		      tiaraDimensions.distTargetEndA) / 2);
  G4double radialHalf(0.5 * (tiaraDimensions.widthExperiment/2 - 
			     tiaraDimensions.pipeRadius));
  G4double lateralHalf(0.5 * (tiaraDimensions.pipeRadius + 
			      tiaraDimensions.widthExperiment/2));
  G4double radialShift(tiaraDimensions.pipeRadius + radialHalf);
  G4double lateralShift(tiaraDimensions.widthExperiment/2 - lateralHalf);

  G4Box *ironBbox = new G4Box("ironBbox",
			      lateralHalf,
			      radialHalf,
			      halfHight);
  
  G4LogicalVolume *logIronBbox = 
    new G4LogicalVolume(ironBbox, 
			mfac.GetMaterial("iron"),
			"IronBbox");
  
  G4VisAttributes* pIronAVis = new 
    G4VisAttributes(G4Colour(0.1, 0.1, 1.0));
  //  pIronAVis->SetForceSolid(true);
  logIronBbox->SetVisAttributes(pIronAVis);
  
  TiaraPart tiaraPart;
  tiaraPart.rot = 0;
  tiaraPart.logVol = logIronBbox;
  tiaraPart.pos.setZ(tiaraDimensions.targetPosZ +
		     tiaraDimensions.distTargetEndA + halfHight);

  // upper shield
  tiaraPart.pos.setX(lateralShift);
  tiaraPart.pos.setY(radialShift);
  fTiaraParts.push_back(tiaraPart);

  // lower shield
  tiaraPart.pos.setX(-1*lateralShift);
  tiaraPart.pos.setY(-1*radialShift);
  fTiaraParts.push_back(tiaraPart);

  // right shield
  tiaraPart.rot = new G4RotationMatrix();
  tiaraPart.rot->set(G4ThreeVector(0, 0, 1), 90 * deg);  
  tiaraPart.pos.setX(-1* radialShift);
  tiaraPart.pos.setY(lateralShift);
  fTiaraParts.push_back(tiaraPart);

  // left shield
  //  tiaraPart.rot->set(G4ThreeVector(0, 0, 1), -90);  
  tiaraPart.pos.setX(radialShift);
  tiaraPart.pos.setY(-1*lateralShift);
  fTiaraParts.push_back(tiaraPart);
  

}
*/

TiaraIronShieldB::
TiaraIronShieldB(TiaraMaterials &mfac,
		 const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetExperiment - 
		      tiaraDimensions.distTargetEndA) / 2);
  G4double halfWidth(tiaraDimensions.widthExperiment/2); 

  G4Box *ironBbox = new G4Box("ironBbox",
			      halfWidth,
			      halfWidth,
			      halfHight);
  
  G4LogicalVolume *logIronBbox = 
    new G4LogicalVolume(ironBbox, 
			mfac.GetMaterial("iron"),
			"IronBbox");
  
  G4VisAttributes* pIronAVis = new 
    G4VisAttributes(G4Colour(0.1, 0.1, 1.0));
  //  pIronAVis->SetForceSolid(true);
  logIronBbox->SetVisAttributes(pIronAVis);

  G4Tubs *tub = new G4Tubs("holeIronB",
			   0,
			   tiaraDimensions.pipeRadius,
			   halfHight,
			   0,
			   360 * deg);
  
  G4LogicalVolume *logHole = new G4LogicalVolume(tub, 
						 mfac.GetMaterial("air"),
						 "logHoleIronB");
  
  new G4PVPlacement(0,  
		    G4ThreeVector(0, 0, 0),
		    logHole,
		    "holeIronB", 
		    logIronBbox, 
		    false, 
		    0);
  
  TiaraPart tiaraPart;
  tiaraPart.rot = 0;
  tiaraPart.logVol = logIronBbox;
  tiaraPart.pos.setZ(tiaraDimensions.targetPosZ +
		     tiaraDimensions.distTargetEndA + halfHight);
  tiaraPart.pos.setX(0);
  tiaraPart.pos.setY(0);
  fTiaraParts.push_back(tiaraPart);

}

TiaraIronShieldB::~TiaraIronShieldB(){
}

TiaraParts TiaraIronShieldB::GetParts(){
  return fTiaraParts;
}
