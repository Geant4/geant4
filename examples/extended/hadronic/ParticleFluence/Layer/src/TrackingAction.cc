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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "Run.hh"

const std::array< G4String, TrackingAction::fkNumberScoringVolumes >
TrackingAction::fkArrayScoringVolumeNames = { "layer" };

const std::array< G4String, TrackingAction::fkNumberKinematicRegions >
TrackingAction::fkArrayKinematicRegionNames = { "", "below 20 MeV", "above 20 MeV" };

const std::array< G4String, TrackingAction::fkNumberParticleTypes >
TrackingAction::fkArrayParticleTypeNames = { "all", "electron", "gamma", "muon", "neutrino",
                                             "pion", "neutron", "proton", "ion", "otherMeson",
                                             "otherBaryon" };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TrackingAction::GetIndex( const G4int iScoringVolume, const G4int iKinematicRegion,
                                const G4int iParticleType ) {
  G4int index = -1;
  if ( iScoringVolume >= 0    &&  iScoringVolume < fkNumberScoringVolumes      &&
       iKinematicRegion >= 0  &&  iKinematicRegion < fkNumberKinematicRegions  &&
       iParticleType >= 0     &&  iParticleType < fkNumberParticleTypes           ) {
    index = iScoringVolume * fkNumberKinematicRegions * fkNumberParticleTypes +
                                     iKinematicRegion * fkNumberParticleTypes +
                                                               iParticleType;
  }
  return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction() : G4UserTrackingAction() {
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::Initialize() {
  // Initialization needed at the beginning of each Run  
  fArrayMultiplicities.fill( 0 );
  fArraySumKineticEnergies.fill( 0.0 );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction( const G4Track* aTrack ) {
  // This method is called not only once when a particle is created, 
  // but also each time it is resumed, in the case the track gets suspended,
  // as it happens in the case of neutrons with _HP Physics Lists.
  // To be sure that we collect information about a track one and only once,
  // we require that the current step be the first one.
  if ( aTrack == nullptr  ||
       aTrack->GetCurrentStepNumber() != 0  ||
       aTrack->GetDefinition() == nullptr  ||
       aTrack->GetLogicalVolumeAtVertex() == nullptr  ||
       aTrack->GetLogicalVolumeAtVertex()->GetName() != "logicLayer" ) return;
  G4int iScoringVolume = 0;
  // Three kinematical regions:  [0] : any value ;  [1] : below 20 MeV ;  [2] : above 20 MeV
  G4int iKinematicRegion = aTrack->GetKineticEnergy() < 20.0 ? 1 : 2;
  G4int absPdg = std::abs( aTrack->GetDefinition()->GetPDGEncoding() );  
  G4int iParticleType = -1;
  if      ( absPdg == 11 ) iParticleType = 1;    // electron (and positron)
  else if ( absPdg == 22 ) iParticleType = 2;    // gamma
  else if ( absPdg == 13 ) iParticleType = 3;    // muons (mu- and mu+)
  else if ( absPdg == 12 || absPdg == 14 || absPdg == 16 ) iParticleType = 4;
    // neutrinos (and anti-neutrinos), all flavors
  else if ( absPdg == 111 || absPdg == 211 ) iParticleType = 5;  // (charged) pions
  else if ( absPdg == 2112 ) iParticleType = 6;  // neutron (and anti-neutron)
  else if ( absPdg == 2212 ) iParticleType = 7;  // proton  (and anti-proton)
  else if ( G4IonTable::IsIon( aTrack->GetDefinition() ) ||
            G4IonTable::IsAntiIon( aTrack->GetDefinition() ) ) iParticleType = 8;
    // ions (and anti-ions)
  else if ( absPdg < 1000 ) iParticleType = 9;   // other mesons (e.g. kaons)
    // (Note: this works in most cases, but not always!)
  else if ( absPdg > 1000 ) iParticleType = 10;  // other baryons (e.g. hyperons,
                                                 // anti-hyperons, etc.)
  // Consider the specific case : scoring volume, kinematic region and particle type
  G4int index = GetIndex( iScoringVolume, iKinematicRegion, iParticleType );
  ++fArrayMultiplicities[index];
  fArraySumKineticEnergies[index] += aTrack->GetKineticEnergy();
  // Consider the "all" particle case, with the same scoring volume and kinematic region
  index = GetIndex( iScoringVolume, iKinematicRegion, 0 );
  ++fArrayMultiplicities[index];
  fArraySumKineticEnergies[index] += aTrack->GetKineticEnergy();
  // Consider the "any" kinematic region case, with the same scoring volume and particle type
  index = GetIndex( iScoringVolume, 0, iParticleType );
  ++fArrayMultiplicities[index];
  fArraySumKineticEnergies[index] += aTrack->GetKineticEnergy();
  // Consider the "any" kinematic region and "all" particle, with the same scoring volume
  index = GetIndex( iScoringVolume, 0, 0 );
  ++fArrayMultiplicities[index];
  fArraySumKineticEnergies[index] += aTrack->GetKineticEnergy();
  if ( fRunPtr ) {
    fRunPtr->SetTrackingArray1( fArrayMultiplicities );
    fRunPtr->SetTrackingArray2( fArraySumKineticEnergies );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction( const G4Track* /* aTrack  */ ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
