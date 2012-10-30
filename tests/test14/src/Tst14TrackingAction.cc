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
#include "Tst14TrackingAction.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"

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
    G4ThreeVector newPol(std::cos(phi),std::sin(phi),0.); 
    G4ParticleMomentum aParticleDirection = aParticle->GetMomentumDirection(); 
    newPol.rotateUz(aParticleDirection);
    ((G4Track*)aTrack)->SetPolarization(newPol);
  }
}
