///////////////////////
//G4RTTrackingAction.cc
///////////////////////


#include "G4RTTrackingAction.hh"
#include "G4RayTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4ios.hh"


void G4RTTrackingAction :: PreUserTrackingAction(const G4Track* aTrack)
{
  G4RayTrajectory* aTrajectory=new G4RayTrajectory;
  fpTrackingManager->SetTrajectory(aTrajectory);

}







