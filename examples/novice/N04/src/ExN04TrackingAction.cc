

#include "ExN04TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

void ExN04TrackingAction::PreUserTrackingAction()
{
  G4TrackingManager* trackingManager = GetOmnipotentTrackingManager();
  G4Track* aTrack = trackingManager->GetTrack();

  // Create trajectory only for primaries
  if(aTrack->GetParentID()==0)
  { trackingManager->SetStoreTrajectory(true); }
  else
  { trackingManager->SetStoreTrajectory(false); }
}


