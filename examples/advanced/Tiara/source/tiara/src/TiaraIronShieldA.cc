#include "TiaraIronShieldA.hh"
#include "TiaraDimensions.hh"
#include "TiaraMaterials.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "TiaraPart.hh"
#include "G4VisAttributes.hh"


TiaraIronShieldA::
TiaraIronShieldA(TiaraMaterials &mfac,
		 const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetEndA - 
		      tiaraDimensions.distTargetWall) / 2);
  G4Tubs *ironAtub = new G4Tubs("ironAtub",
				   tiaraDimensions.pipeRadius,
				   tiaraDimensions.radiusIronA,
				   halfHight,
				   0,
				   360 * deg);
  G4LogicalVolume *logIronA = 
    new G4LogicalVolume(ironAtub, 
			mfac.GetMaterial("iron"),
			"IronA");

  G4VisAttributes* pIronAVis = new 
    G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  //  pIronAVis->SetForceSolid(true);
  logIronA->SetVisAttributes(pIronAVis);
  
  TiaraPart tiaraPart;
  tiaraPart.logVol = logIronA;
  tiaraPart.pos.setX(0);
  tiaraPart.pos.setY(0);
  tiaraPart.pos.setZ(tiaraDimensions.targetPosZ +  
		     tiaraDimensions.distTargetWall + halfHight);
  fTiaraParts.push_back(tiaraPart);
}

TiaraIronShieldA::~TiaraIronShieldA(){
}

TiaraParts TiaraIronShieldA::GetParts(){
  return fTiaraParts;
}
