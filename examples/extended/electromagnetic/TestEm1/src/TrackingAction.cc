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
// $Id: TrackingAction.cc,v 1.8 2004/07/23 15:39:40 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(PrimaryGeneratorAction* prim, RunAction* run, 
                               HistoManager* histo)
:primary(prim), runAction(run), histoManager(histo)
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
  G4double Trleng = aTrack->GetTrackLength();
    
  if (aTrack->GetDefinition()->GetPDGCharge() == 0.) {
     runAction->CountTraks0(1); 
     runAction->CountSteps0(nbSteps);
  
  } else {
     runAction->CountTraks1(1); 
     runAction->CountSteps1(nbSteps);
     histoManager->FillHisto(1,Trleng);
     histoManager->FillHisto(2,(float)nbSteps);
  }
  
  //csda and projected ranges for primary particle
  if (aTrack->GetTrackID() == 1) {
    runAction->AddCsdaRange(Trleng);
    G4ThreeVector vertex = primary->GetParticleGun()->GetParticlePosition();    
    G4ThreeVector position = aTrack->GetPosition() - vertex;      
    runAction->AddProjRange(position.x());
    runAction->AddTransvDev(position.y());
    runAction->AddTransvDev(position.z());
  }        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

