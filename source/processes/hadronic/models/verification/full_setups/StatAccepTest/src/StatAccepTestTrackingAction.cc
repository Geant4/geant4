#include "StatAccepTestTrackingAction.hh"
#include "StatAccepTestAnalysis.hh"


StatAccepTestTrackingAction::StatAccepTestTrackingAction() {}


StatAccepTestTrackingAction::~StatAccepTestTrackingAction() {}


void StatAccepTestTrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {
  StatAccepTestAnalysis::getInstance()->infoTrack( aTrack );
}


void StatAccepTestTrackingAction::PostUserTrackingAction( const G4Track* aTrack ) {}

