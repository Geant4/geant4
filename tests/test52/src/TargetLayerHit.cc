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
