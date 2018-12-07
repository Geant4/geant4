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
//
/// \file Par02FastSimModelEMCal.cc
/// \brief Implementation of the Par02FastSimModelEMCal class

#include "Par02FastSimModelEMCal.hh"
#include "Par02EventInformation.hh"
#include "Par02PrimaryParticleInformation.hh"
#include "Par02Smearer.hh"
#include "Par02Output.hh"

#include "G4Track.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "g4root.hh"

#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelEMCal::Par02FastSimModelEMCal( G4String aModelName, 
  G4Region* aEnvelope, Par02DetectorParametrisation::Parametrisation aType ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( aType ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelEMCal::Par02FastSimModelEMCal( G4String aModelName, 
                                                G4Region* aEnvelope ) : 
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( Par02DetectorParametrisation::eCMS ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelEMCal::Par02FastSimModelEMCal( G4String aModelName ) :
  G4VFastSimulationModel( aModelName ), fCalculateParametrisation(), 
  fParametrisation( Par02DetectorParametrisation::eCMS ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelEMCal::~Par02FastSimModelEMCal() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par02FastSimModelEMCal::IsApplicable( 
  const G4ParticleDefinition& aParticleType ) {
  // Applicable for electrons, positrons, and gammas
  return &aParticleType == G4Electron::Definition()  ||
         &aParticleType == G4Positron::Definition()  ||
         &aParticleType == G4Gamma::Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par02FastSimModelEMCal::ModelTrigger( const G4FastTrack& /*aFastTrack*/ ) {
  return true;  // No kinematical restrictions to apply the parametrisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02FastSimModelEMCal::DoIt( const G4FastTrack& aFastTrack,
                                   G4FastStep& aFastStep ) {
  //G4cout << " ________EMCal model triggered _________" << G4endl;

  // Kill the parameterised particle at the entrance of the electromagnetic calorimeter
  aFastStep.KillPrimaryTrack();
  aFastStep.ProposePrimaryTrackPathLength( 0.0 );
  G4double Edep = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();

  // Consider only primary tracks (do nothing else for secondary e-, e+, gammas)
  G4ThreeVector Pos = aFastTrack.GetPrimaryTrack()->GetPosition();
  if ( ! aFastTrack.GetPrimaryTrack()->GetParentID() ) {
    Par02EventInformation* info = (Par02EventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();
    if ( info->GetDoSmearing() ) {
      // Smearing according to the electromagnetic calorimeter resolution
      G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
      G4double res = fCalculateParametrisation->GetResolution( 
        Par02DetectorParametrisation::eEMCAL, fParametrisation, Porg.mag() );
      G4double eff = fCalculateParametrisation->GetEfficiency( 
        Par02DetectorParametrisation::eEMCAL, fParametrisation, Porg.mag() );
      G4double Esm;
      Esm = std::abs( Par02Smearer::Instance()->
                        SmearEnergy( aFastTrack.GetPrimaryTrack(), res ) );
      Par02Output::Instance()->FillHistogram( 1, (Esm/MeV) / (Edep/MeV) );
      // Setting the values of Pos, Esm, res and eff
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetEMCalPosition( Pos );
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetEMCalEnergy( Esm );
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetEMCalResolution( res );
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetEMCalEfficiency( eff );
      // The (smeared) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the electromagnetic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( Esm );
    } else {
      // No smearing: simply setting the value of Edep
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetEMCalEnergy( Edep );
      // The (initial) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the electromagnetic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( Edep );
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

