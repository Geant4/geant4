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
#include "Tst14TrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "globals.hh"

void Tst14TrackingAction::PreUserTrackingAction(const G4Track* aTrack )
{
  //  G4TrackingManager* trackingManager = fpTrackingManager;
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
