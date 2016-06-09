//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: TrackingAction.cc,v 1.5 2003/10/28 10:27:26 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
// 
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "RunAction.hh"

#include "G4Track.hh"
 
#ifdef G4ANALYSIS_USE
 #include "AIDA/IHistogram1D.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(RunAction* RunAct)
:runAction(RunAct)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track*)
{
  //  G4cout << "ID= " << aTrack->GetTrackID() << "  e(MeV)= " 
  //         << aTrack->GetDynamicParticle()->GetKineticEnergy()/MeV << "  "
  //         << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName()
  //         << G4endl;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //increase nb of processed tracks 
  //count nb of steps of this track
  G4int   nbSteps = aTrack->GetCurrentStepNumber();

  if (aTrack->GetDefinition()->GetPDGCharge() == 0.) {
     runAction->CountTraks0(1); 
     runAction->CountSteps0(nbSteps);
  
  } else {

     runAction->CountTraks1(1); 
     runAction->CountSteps1(nbSteps);

#ifdef G4ANALYSIS_USE
     G4double Trleng = aTrack->GetTrackLength();
     runAction->GetHisto(0)->fill(Trleng);
     runAction->GetHisto(1)->fill((float)nbSteps);
#endif	
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

