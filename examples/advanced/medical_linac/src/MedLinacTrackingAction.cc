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


#include "MedLinacTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
  
#include "MedLinacDetectorConstruction.hh" 
#include "MedLinacTargetAndFilterDecorator.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

void MedLinacTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
 {

  
    // global geometry navigator
    gNavigator = G4TransportationManager::GetTransportationManager()
      ->GetNavigatorForTracking();


   // Create trajectory only for primaries

      // check if particle is gamma
   G4ParticleDefinition* particleType = aTrack->GetDefinition();
   G4ThreeVector direction = aTrack->GetMomentumDirection();
   G4double theta  = std::acos(std::abs(direction.z()));

   //if(aTrack->GetParentID()==0)
   //fpTrackingManager->SetStoreTrajectory(true); 

       if(particleType==G4Gamma::GammaDefinition()) {
	if(theta>=20.*deg){
        // check if particle is in target
        G4ThreeVector pos = aTrack->GetPosition();
        G4ThreeVector *ptr = NULL;
        G4VPhysicalVolume *theVolume;
	theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,false);
        if(theVolume->GetName() == "targetA"||theVolume->GetName() == "targetB") {
    
	  fpTrackingManager->SetStoreTrajectory(false); }
	}
}
       else
	 {fpTrackingManager->SetStoreTrajectory(true);} 
}
