#define hTestTrackingAction_CPP

//---------------------------------------------------------------------------
//
// ClassName:   hTestTrackingAction
//  
// Description: Implementation file for MC truth.
//
// Author:      V.Ivanchenko 17/03/01
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// HARP includes
#include "Simulation/hTestTrackingAction.h"
#include "HarpEvent/MonteCarlo/McParticle.h"
#include "HarpEvent/MonteCarlo/McVertex.h"

// G4 includes
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestTrackingAction::hTestTrackingAction(hTestRunAction* run):
  theRun(run)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestTrackingAction::~hTestTrackingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  if(1 < theRun->GetVerbose()) {
    G4cout << "hTestTrackingAction: Next track #" 
           << aTrack->GetTrackID() << G4endl;
  }

  G4bool primary = false;
  if(0 == aTrack->GetParentID()) primary = true;

  //Save primary parameters

  if(primary) {

    G4double kinE = aTrack->GetKineticEnergy();
    theRun->SaveToTuple("TKIN", kinE/MeV);      

    G4ParticleDefinition* particle = aTrack->GetDefinition();
    G4DynamicParticle* dynParticle = aTrack->GetDynamicParticle();

    G4double mass = 0.0;
    if(particle) {
      mass = particle->GetPDGMass();
      theRun->SaveToTuple("MASS", mass/MeV);      
      theRun->SaveToTuple("CHAR",(particle->GetPDGCharge())/eplus);      
    }

    G4ThreeVector pos = aTrack->GetPosition();
    theRun->SaveToTuple("X0",(pos.x())/mm);
    theRun->SaveToTuple("Y0",(pos.y())/mm);
    theRun->SaveToTuple("Z0",(pos.z())/mm);
  

    if(1 < theRun->GetVerbose()) {
      G4cout << "hTestTrackingAction: kinE(MeV)= " << kinE/MeV
           << "; m(MeV)= " << mass/MeV   
           << "; pos= " << pos << ";  dir= " << dir << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....








