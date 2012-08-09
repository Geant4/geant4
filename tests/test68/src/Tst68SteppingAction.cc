#include "Tst68SteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"


Tst68SteppingAction::Tst68SteppingAction() : 
  totalEdepAllParticles( 0.0 ),
  primaryParticleId( 0 ), primaryParticleEnergy( 0.0 ),
  isFirstStepOfTheEvent( true ) {}


Tst68SteppingAction::~Tst68SteppingAction() {}


void Tst68SteppingAction::UserSteppingAction( const G4Step * theStep ) {

  // Store the information about the ID and the kinetic energy 
  // of the primary particle, at the first step of the first event.
  // Notice that for the kinetic energy, we are considering the 
  // "pre-step" point of such first step.
  // Be careful that the condition "primaryParticleId == 0"
  // cannot be used to check if it is the first step, because
  // ions have PDG code = 0.
  if ( isFirstStepOfTheEvent ) {
    if ( theStep->GetTrack()->GetParentID() == 0 ) {
      primaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      primaryParticleEnergy = theStep->GetPreStepPoint()->GetKineticEnergy();
      isFirstStepOfTheEvent = false;
    }
  }

  // Consider the energy depositions only inside the calorimeter, i.e.
  // everywhere but in the "expHall" volume.
  if ( theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() != "expHall" ) {

    G4double edep = theStep->GetTotalEnergyDeposit() * theStep->GetTrack()->GetWeight();
    // Multiply the energy deposit with the weight of the track,
    // to allow the use of biasing.

    if ( edep > 0.001*eV ) totalEdepAllParticles += edep;
  }

}


void Tst68SteppingAction::reset() {
  totalEdepAllParticles = 0.0;
  primaryParticleId = 0;
  primaryParticleEnergy = 0.0;
  isFirstStepOfTheEvent = true;
}

