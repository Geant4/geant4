#include "TiaraConcreteShieldA.hh"
#include "TiaraDimensions.hh"
#include "TiaraMaterials.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"


TiaraConcreteShieldA::
TiaraConcreteShieldA(TiaraMaterials &mfac,
		     const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetEndA - 
		      tiaraDimensions.distTargetWall) / 2);
  G4Tubs *concretAtub = new G4Tubs("concretAtub",
				   tiaraDimensions.radiusIronA,
				   tiaraDimensions.worldHalfWidth,
				   halfHight,
				   0,
				   360 * deg);
  G4LogicalVolume *logConcreteA = 
    new G4LogicalVolume(concretAtub, 
			mfac.GetMaterial("concrete"),
			"ConcreteA");

  G4VisAttributes* pConcreteAVis = new 
    G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  //  pConcreteAVis->SetForceSolid(true);
  logConcreteA->SetVisAttributes(pConcreteAVis);



  TiaraPart tiaraPart;
  tiaraPart.rot = 0;
  tiaraPart.logVol = logConcreteA;
  tiaraPart.pos.setX(0);
  tiaraPart.pos.setY(0);
  tiaraPart.pos.setZ(tiaraDimensions.targetPosZ +  
		     tiaraDimensions.distTargetWall + halfHight);
  fTiaraParts.push_back(tiaraPart);
}

TiaraConcreteShieldA::~TiaraConcreteShieldA(){
}

TiaraParts TiaraConcreteShieldA::GetParts(){
  return fTiaraParts;
}
