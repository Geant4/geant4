
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2TrackingAction.cc,v 1.2 1999-12-15 14:49:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2TrackingAction.hh"
#include "Em2RunAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2TrackingAction::Em2TrackingAction(Em2RunAction* run)
:Em2Run(run)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //count total track length
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();
  G4double TrLeng = aTrack->GetTrackLength();
  
  Em2Run->fillPerTrack(charge,TrLeng);     
}


