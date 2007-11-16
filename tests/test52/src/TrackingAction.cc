
#include "TrackingAction.hh"
#include "TrackInformation.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"


G4Allocator<TrackInformation> TrackInformationAllocator;


void TrackingAction::PreUserTrackingAction(const G4Track*) {

  fpTrackingManager -> SetUserTrackInformation(new TrackInformation);
}


void TrackingAction::PostUserTrackingAction(const G4Track*) {
 
}
