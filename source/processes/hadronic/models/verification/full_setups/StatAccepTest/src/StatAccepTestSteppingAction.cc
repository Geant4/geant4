#include "StatAccepTestSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"


StatAccepTestSteppingAction::StatAccepTestSteppingAction() : 
  totalEdepAllParticles( 0.0 ),
  primaryParticleId( 0 ), primaryParticleEnergy( 0.0 ) {}


StatAccepTestSteppingAction::~StatAccepTestSteppingAction() {}


void StatAccepTestSteppingAction::UserSteppingAction(const G4Step * theStep) {

  if ( primaryParticleId == 0 ) {
    if ( theStep->GetTrack()->GetParentID() == 0 ) {
      // PDG ID: e- = 11, mu- = 13, pion- = -211, pion+ = 211, p = 2212, n = 2112.
      primaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      primaryParticleEnergy = theStep->GetTrack()->GetKineticEnergy();
    }
  }

  G4double edep = theStep->GetTotalEnergyDeposit();
  if ( edep < 0.001*eV ) return;
 
  totalEdepAllParticles += edep;

}


void StatAccepTestSteppingAction::reset() {
  totalEdepAllParticles = 0.0;
  primaryParticleId = 0;
  primaryParticleEnergy = 0.0;
}

