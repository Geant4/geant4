#include "ML2TrackingAction.h"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
  
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"


void CML2TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
 {

    // global geometry navigator
    //gNavigator = G4TransportationManager::GetTransportationManager()
    //  ->GetNavigatorForTracking();

   // Create trajectory only for primaries

      // check if particle is gamma
   //G4ParticleDefinition* particleType = aTrack->GetDefinition();
   //G4ThreeVector direction = aTrack->GetMomentumDirection();
   //G4double theta  = std::acos(std::abs(direction.z()));
   if(aTrack->GetParentID()!=0)
   {
	   this->fpTrackingManager->SetStoreTrajectory(true); 
   }

//       if(particleType==G4Gamma::GammaDefinition()) {
//	if(theta>=20.*deg){
//        // check if particle is in target
//        G4ThreeVector pos = aTrack->GetPosition();
//        G4ThreeVector *ptr = NULL;
//        G4VPhysicalVolume *theVolume;
//	theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,false);
//        if(theVolume->GetName() == "targetA"||theVolume->GetName() == "targetB") {
//    
//	  fpTrackingManager->SetStoreTrajectory(false); }
//	}
//}
//       else
//	 {fpTrackingManager->SetStoreTrajectory(true);} 
}


CML2TrackingAction::~CML2TrackingAction(void)
{
}

CML2TrackingAction::CML2TrackingAction(void)
{
}

