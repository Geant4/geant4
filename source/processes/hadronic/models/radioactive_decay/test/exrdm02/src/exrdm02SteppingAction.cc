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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4ios.hh"

//#include "CLHEP/Hist/TupleManager.h"
//#include "CLHEP/Hist/HBookFile.h"
//#include "CLHEP/Hist/Histogram.h"
//#include "CLHEP/Hist/Tuple.h"

#include "exrdm02SteppingAction.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"

#include <vector>

extern G4bool drawEvent;

extern std::vector<G4String> Particles;
extern std::vector<G4double> Energies;
extern std::vector<G4double> Weights;
extern std::vector<G4double> Times;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdm02SteppingAction::exrdm02SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdm02SteppingAction::~exrdm02SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdm02SteppingAction::UserSteppingAction(const G4Step* fStep) 
{
  const G4SteppingManager* pSM = fpSteppingManager;
  G4Track* fTrack = pSM->GetTrack();
  //G4Step* fStep = pSM->GetStep();
  //G4int TrackID = fTrack->GetTrackID();
  G4int StepNo = fTrack->GetCurrentStepNumber();
  if(StepNo >= 10000) fTrack->SetTrackStatus(fStopAndKill);
  
  //  cout << fTrack->GetGlobalTime()/s <<"  "<< fTrack->GetLocalTime()/s << " " <<fTrack->GetProperTime()/s << endl;
  //cout << fStep->GetPreStepPoint()->GetGlobalTime() /s << endl; 

  if (StepNo == 1) {
    Particles.push_back ( fTrack->GetDefinition()->GetParticleName() );
    Energies.push_back ( fStep->GetPreStepPoint()->GetKineticEnergy()/keV );
    Weights.push_back ( fStep->GetPreStepPoint()->GetWeight() );
    Times.push_back((fStep->GetPreStepPoint()->GetGlobalTime() - fStep->GetPreStepPoint()->GetLocalTime()) / s );
    drawEvent = true;
    
  }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




