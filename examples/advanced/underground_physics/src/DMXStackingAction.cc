//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// StackingAction program
// --------------------------------------------------------------

#include "DMXStackingAction.hh"

#include "DMXStackingActionMessenger.hh"
#include "DMXDetectorConstruction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"



DMXStackingAction::DMXStackingAction() {

  // new messenger
  theMessenger = new DMXStackingActionMessenger(this);

  // messenger defaults
  killGammasFlag  = 0;

  // global geometry navigator
  gNavigator = G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking();
}


DMXStackingAction::~DMXStackingAction() {
  
  delete theMessenger; 
}


G4ClassificationOfNewTrack DMXStackingAction::ClassifyNewTrack 
(const G4Track* aTrack) {


  static G4int gammasKilled = 0;

  G4ClassificationOfNewTrack classification = fWaiting;

  // Kill gammas from neutrons in concrete wall
  if(killGammasFlag) {
    // check if particle is gamma
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    if(particleType==G4Gamma::GammaDefinition()) {
      // check if particle is in world_phys
      G4ThreeVector pos = aTrack->GetPosition();
      G4ThreeVector *ptr = NULL;
      G4VPhysicalVolume *theVolume;
      theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,false);
      if(theVolume->GetName() == "world_phys") {
	classification = fKill;
	G4cout << "     Gammas killed in concrete wall: "
	       << ++gammasKilled << G4endl;
      }
    }
  }

  return classification;

}


void DMXStackingAction::NewStage() {;}

    
void DMXStackingAction::PrepareNewEvent() {;}




