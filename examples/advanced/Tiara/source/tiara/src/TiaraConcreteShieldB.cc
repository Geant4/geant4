// $Id: TiaraConcreteShieldB.cc,v 1.2 2003-06-16 17:06:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "TiaraConcreteShieldB.hh"
#include "TiaraDimensions.hh"
#include "TiaraMaterials.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"


TiaraConcreteShieldB::
TiaraConcreteShieldB(TiaraMaterials &mfac,
		     const TiaraDimensions &tiaraDimensions){
  G4double halfHight((tiaraDimensions.distTargetEndB - 
		      tiaraDimensions.distTargetEndA) / 2);
  G4double radialHalf(0.5 * (tiaraDimensions.worldHalfWidth - 
			     tiaraDimensions.widthExperiment/2));
  G4double lateralHalf(0.5 * (tiaraDimensions.widthExperiment/2 + 
			      tiaraDimensions.worldHalfWidth));
  G4double radialShift(tiaraDimensions.widthExperiment/2 + radialHalf);
  G4double lateralShift(tiaraDimensions.worldHalfWidth - lateralHalf);

  G4Box *concretBbox = new G4Box("concretBbox",
				 lateralHalf,
				 radialHalf,
				 halfHight);
  
  G4LogicalVolume *logConcreteBbox = 
    new G4LogicalVolume(concretBbox, 
			mfac.GetMaterial("concrete"),
			"ConcreteBbox");
  
  G4VisAttributes* pConcreteAVis = new 
    G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  //  pConcreteAVis->SetForceSolid(true);
  logConcreteBbox->SetVisAttributes(pConcreteAVis);
  
  TiaraPart tiaraPart;
  tiaraPart.rot = 0;
  tiaraPart.logVol = logConcreteBbox;
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

TiaraConcreteShieldB::~TiaraConcreteShieldB(){
}

TiaraParts TiaraConcreteShieldB::GetParts(){
  return fTiaraParts;
}
