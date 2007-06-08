#include "DetectorSliceTrackingAction.hh"
#include "DetectorSliceAnalysis.hh"
#include "G4Track.hh"


DetectorSliceTrackingAction::DetectorSliceTrackingAction() {}


DetectorSliceTrackingAction::~DetectorSliceTrackingAction() {}


void DetectorSliceTrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {
  // This method is called not only once when a particle is created, 
  // but also each time it is resumed, in the case the track gets
  // suspended, as it happens in the case of neutrons with _HP
  // Physics Lists. To be sure that we collect information about
  // a track one and only once, we require that the current step
  // be the first one.
  if ( aTrack->GetCurrentStepNumber() == 0 ) {
    DetectorSliceAnalysis::getInstance()->infoTrack( aTrack );
  } else {
    //G4cout << " aTrack->GetCurrentStepNumber() = " 
    //	     <<  aTrack->GetCurrentStepNumber() << G4endl;  //***DEBUG***
  }
}


void DetectorSliceTrackingAction::PostUserTrackingAction( const G4Track* ) {}

