

#include "Tst09TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "globals.hh"

void Tst09TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  G4TrackingManager* trackingManager = fpTrackingManager;
  const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle(); 

  // Add polarization only for gamma
  if ( aParticle->GetDefinition() == G4Gamma::Gamma() ){

	// check polarization vector
    // if polarization exists, leave as it is
    if (aParticle->GetPolarization().mag() >0.0) return;

    // set gamma polarization
    //     Isotropic distribution of photons with transverse polarization 
    //     with respect to gamma momentum
    G4double  phi = twopi * G4UniformRand();
    G4ThreeVector newPol(cos(phi),sin(phi),0.); 
    G4ParticleMomentum aParticleDirection = aParticle->GetMomentumDirection(); 
    newPol.rotateUz(aParticleDirection);
    ((G4Track*)aTrack)->SetPolarization(newPol);
  }
}


