// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RTTrackingAction.cc,v 1.2 2000-03-09 15:36:38 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//

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







