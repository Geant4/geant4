// $Id: A01TrackingAction.cc,v 1.1 2002-11-13 07:24:08 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

#include "A01Trajectory.hh"

A01TrackingAction::A01TrackingAction()
{;}

A01TrackingAction::~A01TrackingAction()
{;}

void A01TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(new A01Trajectory(aTrack));
}


