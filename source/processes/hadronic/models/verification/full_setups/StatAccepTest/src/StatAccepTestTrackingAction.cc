#include "StatAccepTestTrackingAction.hh"
#include "StatAccepTestAnalysis.hh"
#include "G4Track.hh"


StatAccepTestTrackingAction::StatAccepTestTrackingAction() {}


StatAccepTestTrackingAction::~StatAccepTestTrackingAction() {}


void StatAccepTestTrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {

  // This method is called not only once when a particle is created, 
  // but also each time it is resumed, in the case the track gets
  // suspended, as it happens in the case of neutrons with _HP
  // Physics Lists. To be sure that we collect information about
  // a track one and only once, we require that the current step
  // be the first one.
  if ( aTrack->GetCurrentStepNumber() == 0 ) {
    StatAccepTestAnalysis::getInstance()->infoTrack( aTrack );
  } else {
    //G4cout << " aTrack->GetCurrentStepNumber() = " 
    //	     <<  aTrack->GetCurrentStepNumber() << G4endl;  //***DEBUG***
  }
}


void StatAccepTestTrackingAction::PostUserTrackingAction( const G4Track* aTrack ) {}

