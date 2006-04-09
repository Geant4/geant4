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


void StatAccepTestTrackingAction::PostUserTrackingAction( const G4Track* aTrack ) {
  // This method is called not only once at the end of the life of
  // a track, but also each time it is suspended, as it happens 
  // in the case of neutrons with _HP Physics Lists. 
  // To be sure that we collect information about a track only once
  // when its life comes to the end, we have to require that its 
  // status is "fStopAndKill".
  if ( aTrack->GetTrackStatus() == fStopAndKill ) {
    StatAccepTestAnalysis::getInstance()->infoTrack( aTrack );
    //G4cout << " POST-USER Tracking Action info : " << G4endl   //***DEBUG***
    //       << "\t track id = " << aTrack->GetTrackID() << G4endl
    //	     << "\t particle = " << aTrack->GetDefinition()->GetParticleName() << G4endl
    //	     << "\t track status = " << aTrack->GetTrackStatus() << G4endl
    //	     << "\t track length = " << aTrack->GetTrackLength() / mm << " mm" << G4endl
    //	     << "\t track Ekin = " << aTrack->GetKineticEnergy() << " MeV" << G4endl
    //	     << "\t volume = " << aTrack->GetVolume()->GetName() << G4endl
    //	     << "\t material = " << aTrack->GetMaterial()->GetName() << G4endl;

  }
}

