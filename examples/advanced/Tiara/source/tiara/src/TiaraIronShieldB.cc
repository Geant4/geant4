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
