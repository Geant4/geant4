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
#include "Tst18SteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"

#include "g4std/vector"

extern G4bool drawEvent;

extern G4std::vector<G4String> Particles;
extern G4std::vector<G4double> Energies;
extern G4std::vector<G4double> Weights;
extern G4std::vector<G4double> Times;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18SteppingAction::Tst18SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18SteppingAction::~Tst18SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst18SteppingAction::UserSteppingAction(const G4Step* fStep) 
{
  const G4SteppingManager* pSM = fpSteppingManager;
  G4Track* fTrack = pSM->GetTrack();
  //G4Step* fStep = pSM->GetStep();
  //G4int TrackID = fTrack->GetTrackID();
  G4int StepNo = fTrack->GetCurrentStepNumber();
  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);
  
  if (StepNo == 1) {
    if (Particles.size() < 100) {
      Particles.push_back ( fTrack->GetDefinition()->GetParticleName() );
      Energies.push_back ( fStep->GetPreStepPoint()->GetKineticEnergy()/keV );
      Weights.push_back ( fStep->GetPreStepPoint()->GetWeight() );
      Times.push_back((fStep->GetPreStepPoint()->GetGlobalTime() - fStep->GetPreStepPoint()->GetLocalTime()) / s );
      drawEvent = true;
    }
  }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




