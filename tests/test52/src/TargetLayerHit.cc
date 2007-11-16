
#include "TargetLayerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


TargetLayerHit::TargetLayerHit() :
   volName(""),
   energyDeposit(0),
   positionPre(0),
   positionPost(0),
   particle(""),
   trackID(-1) {

}


TargetLayerHit::TargetLayerHit(const TargetLayerHit& right) :
   G4VHit(right),
   volName(right.volName),
   energyDeposit(right.energyDeposit),
   positionPre(right.positionPre),   
   positionPost(right.positionPost),   
   particle(right.particle),
   trackID(right.trackID),
   secParticles(right.secParticles) {

}


TargetLayerHit::~TargetLayerHit() {

}


const TargetLayerHit& TargetLayerHit::operator= (const TargetLayerHit& right) {

  volName = right.volName;
  energyDeposit = right.energyDeposit;
  positionPre = right.positionPre;
  positionPost = right.positionPost;
  particle = right.particle;
  trackID = right.trackID;
  secParticles = right.secParticles;  

  return *this;
}


int TargetLayerHit::operator==(const TargetLayerHit& right) const {

  return (this == &right) ? 1 : 0;
}


void TargetLayerHit::Draw() {

  G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();

  if(visManager) {

    G4Circle circle(positionPre);
    circle.SetScreenSize(0.01);
    circle.SetFillStyle(G4Circle::filled);

    G4VisAttributes visAttributes(G4Colour::Blue());
    if(trackID > 1) visAttributes.SetColour(G4Colour::Green());
    circle.SetVisAttributes(visAttributes);

    visManager->Draw(circle);
  }
}
