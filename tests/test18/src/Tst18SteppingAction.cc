//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include <vector>

#include "Tst18SteppingAction.hh"
#include "Tst18EventAction.hh"
#include "Tst18RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"


Tst18SteppingAction::Tst18SteppingAction(Tst18RunAction* rA, Tst18EventAction* eA)
 : runAction(rA), eventAction(eA)
{}


Tst18SteppingAction::~Tst18SteppingAction()
{}


void Tst18SteppingAction::UserSteppingAction(const G4Step* fStep) 
{
  const G4SteppingManager* pSM = fpSteppingManager;
  G4Track* fTrack = pSM->GetTrack();
  G4int StepNo = fTrack->GetCurrentStepNumber();
  G4int TrackNo = fTrack->GetTrackID();
  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);

  if (StepNo == 1) {
    if(TrackNo == 1) {
      runAction->FillEventNumber();
    }
    G4StepPoint* prePoint = fStep->GetPreStepPoint();

    runAction->FillParticleName(fTrack->GetDefinition()->GetParticleName() );
    runAction->FillEnergy(prePoint->GetKineticEnergy()/keV);
    runAction->FillWeight(prePoint->GetWeight() ); 
    runAction->FillTime((prePoint->GetGlobalTime() - prePoint->GetLocalTime() )/s);
    eventAction->IncrementParticleNumber();
  } 
}

