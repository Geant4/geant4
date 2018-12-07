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
/// \file Par02FastSimModelHCal.cc
/// \brief Implementation of the Par02FastSimModelHCal class

#include "Par02FastSimModelHCal.hh"
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelHCal::Par02FastSimModelHCal( G4String aModelName, 
  G4Region* aEnvelope, Par02DetectorParametrisation::Parametrisation aType ) :
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( aType ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelHCal::Par02FastSimModelHCal( G4String aModelName, 
                                              G4Region* aEnvelope ) : 
  G4VFastSimulationModel( aModelName, aEnvelope ), fCalculateParametrisation(),
  fParametrisation( Par02DetectorParametrisation::eCMS ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelHCal::Par02FastSimModelHCal( G4String aModelName ) :
  G4VFastSimulationModel( aModelName ), fCalculateParametrisation(), 
  fParametrisation( Par02DetectorParametrisation::eCMS ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02FastSimModelHCal::~Par02FastSimModelHCal() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par02FastSimModelHCal::IsApplicable( const G4ParticleDefinition& aParticleType ) {
  G4bool isOk = false;
  // Applicable to all hadrons, i.e. any particle made of quarks
  if ( aParticleType.GetQuarkContent(1) +
       aParticleType.GetQuarkContent(2) +
       aParticleType.GetQuarkContent(3) +
       aParticleType.GetQuarkContent(4) +
       aParticleType.GetQuarkContent(5) +
       aParticleType.GetQuarkContent(6) != 0 ) {
    isOk = true;
  }
   return isOk;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par02FastSimModelHCal::ModelTrigger( const G4FastTrack& /*aFastTrack*/ ) {
  return true;  // No kinematical restrictions to apply the parametrisation
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par02FastSimModelHCal::DoIt( const G4FastTrack& aFastTrack, 
                                  G4FastStep& aFastStep ) {
  //G4cout << " ________HCal model triggered _________" << G4endl;

  // Kill the parameterised particle at the entrance of the hadronic calorimeter
  aFastStep.KillPrimaryTrack();
  aFastStep.ProposePrimaryTrackPathLength( 0.0 );
  G4double Edep = aFastTrack.GetPrimaryTrack()->GetKineticEnergy();

  // Consider only primary tracks (do nothing else for secondary hadrons)
  G4ThreeVector Pos = aFastTrack.GetPrimaryTrack()->GetPosition();
  if ( ! aFastTrack.GetPrimaryTrack()->GetParentID() ) {
    Par02EventInformation* info = (Par02EventInformation*) 
                            G4EventManager::GetEventManager()->GetUserInformation();
    if ( info->GetDoSmearing() ) {
      // Smearing according to the hadronic calorimeter resolution
      G4ThreeVector Porg = aFastTrack.GetPrimaryTrack()->GetMomentum();
      G4double res = fCalculateParametrisation->
        GetResolution( Par02DetectorParametrisation::eHCAL, 
                       fParametrisation, Porg.mag() );
      G4double eff = fCalculateParametrisation->
        GetEfficiency( Par02DetectorParametrisation::eHCAL, 
                       fParametrisation, Porg.mag() );
      G4double Esm;
      Esm = std::abs( Par02Smearer::Instance()->
                        SmearEnergy( aFastTrack.GetPrimaryTrack(), res ) );
      Par02Output::Instance()->FillHistogram( 2, (Esm/MeV) / (Edep/MeV) );
      // Setting the values of Pos, Esm, res and eff
      Par02PrimaryParticleInformation* primaryInfo= 
         static_cast<Par02PrimaryParticleInformation*>(
           ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) ;
      primaryInfo->SetHCalPosition( Pos );
      primaryInfo->SetHCalEnergy( Esm );
      primaryInfo->SetHCalResolution( res );
      primaryInfo->SetHCalEfficiency( eff );
      // The (smeared) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the hadronic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( Esm );
    } else {
      // No smearing: simply setting the value of Edep
      ( (Par02PrimaryParticleInformation*) ( const_cast< G4PrimaryParticle* >
          ( aFastTrack.GetPrimaryTrack()->GetDynamicParticle()->GetPrimaryParticle() )->
            GetUserInformation() ) )->SetHCalEnergy( Edep );
      // The (initial) energy of the particle is deposited in the step
      // (which corresponds to the entrance of the hadronic calorimeter)
      aFastStep.ProposeTotalEnergyDeposited( Edep );
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

