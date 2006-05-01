#include "StatAccepTestSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "StatAccepTestAnalysis.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"


StatAccepTestSteppingAction::StatAccepTestSteppingAction() : 
  totalEdepAllParticles( 0.0 ),
  primaryParticleId( 0 ), primaryParticleEnergy( 0.0 ) {}


StatAccepTestSteppingAction::~StatAccepTestSteppingAction() {}


void StatAccepTestSteppingAction::UserSteppingAction(const G4Step * theStep) {

  // 10-Apr-2006 : Added energy spectra.
  // We are considering here the kinetic energy spectra of some types
  // of particles when the enter an active layer.
  // If you want instead the spectra for particles exiting, instead of
  // entering, an active layer, then you have to un-comment some lines
  // and comment out some others (see below).
  // Notice that we are not considering here the energy spectra at 
  // production, but at the entrance (or exiting) of active layers.
  // To have more information, we define "negative kinetic energy" when
  // the particle is going backward with respect the beam direction.
  if ( theStep->GetTrack()->GetNextVolume() ) {
    G4TouchableHistory* thePreTouchable = 
      ( G4TouchableHistory* )( theStep->GetPreStepPoint()->GetTouchable() );
    G4TouchableHistory* thePostTouchable = 
      ( G4TouchableHistory* )( theStep->GetPostStepPoint()->GetTouchable() );
    G4String PreVol  = thePreTouchable->GetVolume()->GetName();
    G4String PostVol = thePostTouchable->GetVolume()->GetName();    
    if ( PreVol != "physiActive"  &&  PostVol == "physiActive" ) { // entering active layer
    // if ( PreVol != PostVol  &&  PreVol == "physiActive" ) { // leaving active layer
      G4int replicaNumber = theStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1); // if entering active layer
      // G4int replicaNumber = theStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(1); // if exiting active layer
      G4double kineticEnergy = theStep->GetTrack()->GetKineticEnergy();
      G4double zMomentum = theStep->GetTrack()->GetMomentum().z();
      if ( zMomentum < 0.0 ) {
	// We define "negative kinetic energy" when the particle is 
        // going backward with respect the beam direction.
	kineticEnergy *= -1.0;
      }
      StatAccepTestAnalysis::getInstance()->fillSpectrum( theStep->GetTrack()->GetDefinition(), replicaNumber, kineticEnergy );
    } 
  }

  StatAccepTestAnalysis::getInstance()->infoStep( theStep );
 
  if ( primaryParticleId == 0 ) {
    if ( theStep->GetTrack()->GetParentID() == 0 ) {
      primaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      primaryParticleEnergy = theStep->GetPreStepPoint()->GetKineticEnergy();
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


void StatAccepTestSteppingAction::reset() {
  totalEdepAllParticles = 0.0;
  primaryParticleId = 0;
  primaryParticleEnergy = 0.0;
}

