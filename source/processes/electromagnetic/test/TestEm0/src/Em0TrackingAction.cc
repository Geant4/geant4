
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0TrackingAction.cc,v 1.1 1999-01-08 16:32:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em0TrackingAction.hh"
#include "Em0RunAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0TrackingAction::Em0TrackingAction(Em0RunAction* RunAct)
:runAction(RunAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0TrackingAction::PostUserTrackingAction()
{
  //increase nb of processed tracks 
  //nb of steps of this track
  G4TrackingManager* trackingManager =  GetOmnipotentTrackingManager();
  G4Track* aTrack  = trackingManager -> GetTrack();
  
  G4int nbSteps = aTrack->GetCurrentStepNumber();
  
  if (aTrack->GetDefinition()->GetPDGCharge() == 0.)
       {runAction->CountTraks0(1); runAction->CountSteps0(nbSteps);}
  else {runAction->CountTraks1(1); runAction->CountSteps1(nbSteps);}      
}


