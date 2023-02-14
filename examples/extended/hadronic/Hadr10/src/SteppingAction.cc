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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction() {
  Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::Initialize() {
  // Initialization needed at the beginning of each Run
  fRunPtr = nullptr;  
  fToleranceEPviolations = 1.0*CLHEP::eV;  //***LOOKHERE***
  fPrimaryParticleId = 0;
  fPrimaryParticleInitialKineticEnergy = 0.0;
  fPrimaryParticleInitialTotalEnergy = 0.0;
  fPrimaryParticleInitialMomentum = 0.0;
  fPrimaryParticleInitialBeta = 1.0; 
  fPrimaryParticleInitialGamma = 1.0; 
  fPrimaryParticleInitial3Momentum = G4ThreeVector( 0.0, 0.0, 0.0 );
  fPrimaryParticleInitialPosition = G4ThreeVector( 0.0, 0.0, 0.0 );
  fMaxEkin_deltaMax = 0.0;
  fMaxEtot_deltaMax = 0.0;
  fMaxP_deltaMax = 0.0;
  fMaxPdir_deltaMax = 0.0;
  fMaxMass_deltaMax1 = 0.0;
  fMaxMass_deltaMax2 = 0.0;
  fMaxMass_deltaMax3 = 0.0;
  fMeanMass_deltaMax3 = 0.0;
  fMaxBeta_deltaMax1 = 0.0;
  fMaxBeta_deltaMax2 = 0.0;
  fMaxGamma_deltaMax1 = 0.0;
  fMaxGamma_deltaMax2 = 0.0;
  fMaxGamma_deltaMax3 = 0.0;
  fMaxT_proper_deltaMax = 0.0;
  fMaxT_lab_deltaMax = 0.0;
  fMaxMc_truth_rPos_deltaMax = 0.0;
  fMeanMc_truth_rPos_deltaMax = 0.0;
  fMeanDeltaR_primaryDecay = 0.0;
  fMinDeltaR_primaryDecay = 9999999.9;
  fMaxDeltaR_primaryDecay = -9999999.9;
  fMeanR_primaryDecay = 0.0;
  fMinR_primaryDecay = 9999999.9;
  fMaxR_primaryDecay = -9999999.9;
  fMeanX_primaryDecay = 0.0;
  fMinX_primaryDecay = 9999999.9;
  fMaxX_primaryDecay = -9999999.9;
  fMeanY_primaryDecay = 0.0;
  fMinY_primaryDecay = 9999999.9;
  fMaxY_primaryDecay = -9999999.9;
  fMeanZ_primaryDecay = 0.0;
  fMinZ_primaryDecay = 9999999.9;
  fMaxZ_primaryDecay = -9999999.9;
  fMeanDeltaAngle_primaryDecay = 0.0;
  fMinDeltaAngle_primaryDecay = 9999999.9;
  fMaxDeltaAngle_primaryDecay = -9999999.9;
  fMeanDeltaEkin_primaryDecay = 0.0;
  fMinDeltaEkin_primaryDecay = 9999999.9;
  fMaxDeltaEkin_primaryDecay = -9999999.9;
  fMeanEkin_primaryDecay = 0.0;
  fMinEkin_primaryDecay = 9999999.9;
  fMaxEkin_primaryDecay = -9999999.9;
  fMeanPx_primaryDecay = 0.0;
  fMinPx_primaryDecay = 9999999.9;
  fMaxPx_primaryDecay = -9999999.9;
  fMeanPy_primaryDecay = 0.0;
  fMinPy_primaryDecay = 9999999.9;
  fMaxPy_primaryDecay = -9999999.9;
  fMeanPz_primaryDecay = 0.0;
  fMinPz_primaryDecay = 9999999.9;
  fMaxPz_primaryDecay = -9999999.9;
  fMinUnderestimated_mc_truth_rPos_delta = 9999999.9;
  fMaxOverestimated_mc_truth_rPos_delta = -9999999.9;
  fMeanUnderestimated_mc_truth_rPos_delta = 0.0;
  fMeanOverestimated_mc_truth_rPos_delta = 0.0; 
  fMinUnderestimated_rDeltaPos = 9999999.9;
  fMaxOverestimated_rDeltaPos = -9999999.9;
  fMeanUnderestimated_rDeltaPos = 0.0;
  fMeanOverestimated_rDeltaPos = 0.0;
  fMaxFloat_rDeltaPos_deltaMax = -9999999.9;
  fMeanViolationE_primaryDecay = 0.0;
  fMinViolationE_primaryDecay = 9999999.9;
  fMaxViolationE_primaryDecay = -9999999.9;
  fMeanViolationPx_primaryDecay = 0.0;
  fMinViolationPx_primaryDecay = 9999999.9;
  fMaxViolationPx_primaryDecay = -9999999.9;
  fMeanViolationPy_primaryDecay = 0.0;
  fMinViolationPy_primaryDecay = 9999999.9;
  fMaxViolationPy_primaryDecay = -9999999.9;
  fMeanViolationPz_primaryDecay = 0.0;
  fMinViolationPz_primaryDecay = 9999999.9;
  fMaxViolationPz_primaryDecay = -9999999.9;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction( const G4Step* theStep ) {
  
  // Store the information about the ID and the kinetic energy of the primary particle,
  // at the first step of the first event.
  // Note that for the kinetic energy, we are considering the "pre-step" point of such first step.
  if ( theStep->GetTrack()->GetParentID() == 0  &&
       theStep->GetTrack()->GetCurrentStepNumber() == 1 ) {
    fPrimaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    fPrimaryParticleInitialKineticEnergy = theStep->GetPreStepPoint()->GetKineticEnergy();
    fPrimaryParticleInitialTotalEnergy   = theStep->GetPreStepPoint()->GetTotalEnergy();
    fPrimaryParticleInitial3Momentum     = theStep->GetPreStepPoint()->GetMomentum();
    fPrimaryParticleInitialMomentum      = fPrimaryParticleInitial3Momentum.mag();
    fPrimaryParticleInitialPosition      = theStep->GetPreStepPoint()->GetPosition();
    fPrimaryParticleInitialBeta          = theStep->GetPreStepPoint()->GetBeta();
    fPrimaryParticleInitialGamma         = theStep->GetPreStepPoint()->GetGamma(); 
    // As tolerance for EP violations, consider the max value between the default value
    // and 1 billionth of the initial, primary particle kinetic energy.
    if ( fToleranceEPviolations < fPrimaryParticleInitialKineticEnergy*1.0e-9 ) {
      fToleranceEPviolations = fPrimaryParticleInitialKineticEnergy*1.0e-9;
    }
    // Set the values of this run to the Run object
    if ( fRunPtr ) {
      fRunPtr->SetPrimaryParticleId( fPrimaryParticleId );
      fRunPtr->SetPrimaryParticleInitialKineticEnergy( fPrimaryParticleInitialKineticEnergy );
      fRunPtr->SetPrimaryParticleInitialTotalEnergy( fPrimaryParticleInitialTotalEnergy );
      fRunPtr->SetPrimaryParticleInitialMomentum( fPrimaryParticleInitialMomentum );
      fRunPtr->SetPrimaryParticleInitialBeta( fPrimaryParticleInitialBeta );
      fRunPtr->SetPrimaryParticleInitialGamma( fPrimaryParticleInitialGamma );
      fRunPtr->SetPrimaryParticleInitial3Momentum( fPrimaryParticleInitial3Momentum );
      fRunPtr->SetPrimaryParticleInitialPosition( fPrimaryParticleInitialPosition );
      fRunPtr->SetToleranceEPviolations( ToleranceEPviolations() );
      fRunPtr->SetToleranceDeltaDecayRadius( ToleranceDeltaDecayRadius() );
      fRunPtr->SetIsPreassignedDecayEnabled( IsPreassignedDecayEnabled() );  
      fRunPtr->SetIsBoostToLabEnabled( IsBoostToLabEnabled() );
    }
    // Use the preassigned decay is enabled
    if ( IsPreassignedDecayEnabled()  &&
         ( ! theStep->GetTrack()->GetDefinition()->GetPDGStable() ) ) {
      G4DynamicParticle* dynamicParent =
        const_cast< G4DynamicParticle* >( theStep->GetTrack()->GetDynamicParticle() );
      if ( dynamicParent != nullptr ) {
        G4DecayProducts* decayProducts =
          (G4DecayProducts*)( dynamicParent->GetPreAssignedDecayProducts() );
        if ( decayProducts == nullptr ) {
          G4ParticleDefinition* parentDef = theStep->GetTrack()->GetDefinition();
          G4DecayTable* decayTable =
            ( parentDef == nullptr ? nullptr : parentDef->GetDecayTable() );
          if ( decayTable != nullptr ) {
            G4double parentMass = dynamicParent->GetMass();
            G4VDecayChannel* decayChannel = decayTable->SelectADecayChannel( parentMass );
            if ( decayChannel != nullptr ) {
              decayProducts = decayChannel->DecayIt( parentMass );
              if ( ! decayProducts->IsChecked() ) decayProducts->DumpInfo();
              if ( IsBoostToLabEnabled() ) {
                // boost all decay products to laboratory frame
                decayProducts->Boost( dynamicParent->GetTotalEnergy(),
                                      dynamicParent->GetMomentumDirection() );
              }
            } else {
              decayProducts = new G4DecayProducts( *dynamicParent );
            }
            dynamicParent->SetPreAssignedDecayProducts( decayProducts );
          }
        } else {
          G4cout << "WARNING : already present preassign decay !" << G4endl;
        }
      }
    }
  }

  //G4cout << theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << G4endl;

  // If the primary decays somewhere inside the World volume, get the information about the decay
  if ( theStep->GetTrack()->GetParentID() == 0  &&
       theStep->GetPostStepPoint()->GetProcessDefinedStep() != nullptr  &&
       theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName().find( "Decay" )
       != std::string::npos ) {
    // Get properties of the primary particle when it decays
    
    //--- Get values in different ways and check their consistency ---
    // Kinetic energy of the primary particle at the decay
    const G4double ekin_dynamicParticle =
      theStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy();
    const G4double ekin_track = theStep->GetTrack()->GetKineticEnergy();
    const G4double ekin_postStepPoint = theStep->GetPostStepPoint()->GetKineticEnergy();
    const G4double ekin_deltaMax =
      std::max( std::abs( ekin_dynamicParticle - ekin_track ),
                std::abs( ekin_dynamicParticle - ekin_postStepPoint ) );
    //G4cout << "\t ekin_deltaMax [eV] = " << ekin_deltaMax / CLHEP::eV << G4endl; 
    const G4double ekin_val = ekin_dynamicParticle;  // To be used later
    // Total energy of the primary particle at the decay
    const G4double etot_dynamicParticle =
      theStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy();
    const G4double etot_track = theStep->GetTrack()->GetTotalEnergy();
    const G4double etot_postStepPoint = theStep->GetPostStepPoint()->GetTotalEnergy();
    const G4double etot_deltaMax =
      std::max( std::abs( etot_dynamicParticle - etot_track ),
                std::abs( etot_dynamicParticle - etot_postStepPoint ) );
    //G4cout << "\t etot_deltaMax [eV] = " << etot_deltaMax / CLHEP::eV << G4endl; 
    const G4double etot_val = etot_dynamicParticle;  // To be used later
    // Module of the 3-momentum of the primary particle at the decay
    const G4double p_dynamicParticle =
      theStep->GetTrack()->GetDynamicParticle()->GetMomentum().mag();
    const G4double p_track = theStep->GetTrack()->GetMomentum().mag();
    const G4double p_postStepPoint = theStep->GetPostStepPoint()->GetMomentum().mag();
    const G4double p_deltaMax = std::max( std::abs( p_dynamicParticle - p_track ),
                                          std::abs( p_dynamicParticle - p_postStepPoint ) );
    //G4cout << "\t p_deltaMax [eV] = " << p_deltaMax / CLHEP::eV << G4endl; 
    const G4double p_val = p_dynamicParticle;  // To be used later
    // 3-momentum direction (adimensional) of the primary particle at the decay
    const G4ThreeVector pdir_dynamicParticle =
      theStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection();
    const G4ThreeVector pdir_track = theStep->GetTrack()->GetMomentumDirection();
    const G4ThreeVector pdir_postStepPoint = theStep->GetPostStepPoint()->GetMomentumDirection();
    const G4double pdir_x_deltaMax =
      std::max( std::abs( pdir_dynamicParticle.x() - pdir_track.x() ),
                std::abs( pdir_dynamicParticle.x() - pdir_postStepPoint.x() ) );     
    const G4double pdir_y_deltaMax =
      std::max( std::abs( pdir_dynamicParticle.y() - pdir_track.y() ),
                std::abs( pdir_dynamicParticle.y() - pdir_postStepPoint.y() ) );     
    const G4double pdir_z_deltaMax =
      std::max( std::abs( pdir_dynamicParticle.z() - pdir_track.z() ),
                std::abs( pdir_dynamicParticle.z() - pdir_postStepPoint.z() ) );
    const G4double pdir_deltaMax =
      std::max( std::max( pdir_x_deltaMax, pdir_y_deltaMax ), pdir_z_deltaMax );
    //G4cout << "\t pdir_deltaMax = " << pdir_deltaMax << G4endl; 
    // Mass of the primary particle at the decay
    const G4double mass_dynamicParticle = theStep->GetTrack()->GetDynamicParticle()->GetMass();
    const G4double mass_preStepPoint    = theStep->GetPreStepPoint()->GetMass();
    const G4double mass_postStepPoint   = theStep->GetPostStepPoint()->GetMass();
    const G4double mass_from_etot_ekin  = etot_val - ekin_val;
    const G4double mass_from4mom        = std::sqrt( etot_val*etot_val - p_val*p_val );
    G4double mass_deltaMax1 = std::max( std::abs( mass_dynamicParticle - mass_preStepPoint ),
                                        std::abs( mass_dynamicParticle - mass_postStepPoint ) );
    G4double mass_deltaMax2 = std::abs( mass_dynamicParticle - mass_from_etot_ekin );
    G4double mass_deltaMax3 = std::abs( mass_dynamicParticle - mass_from4mom );
    fMeanMass_deltaMax3 += mass_deltaMax3;
    //G4cout << "\t mass_deltaMax{1,2,3} [eV] = " << mass_deltaMax1 / CLHEP::eV << "\t"
    //       << mass_deltaMax2 / CLHEP::eV << "\t" << mass_deltaMax3 / CLHEP::eV << G4endl;
    const G4double mass_val = mass_dynamicParticle;  // To be used later
    // Lorentz beta of the primary particle at the decay
    // The following line works only for G4 versions >= 10.7
    const G4double beta_dynamicParticle = theStep->GetTrack()->GetDynamicParticle()->GetBeta();
    const G4double beta_postStepPoint   = theStep->GetPostStepPoint()->GetBeta();
    //Before-10.7  const G4double beta_dynamicParticle = beta_postStepPoint;
    const G4double beta_velocity_track = theStep->GetTrack()->GetVelocity() / CLHEP::c_light;
    const G4double beta_velocity_postStepPoint =
      theStep->GetPostStepPoint()->GetVelocity() / CLHEP::c_light;
    const G4double beta_p_over_etot = p_val / etot_val;
    G4double beta_deltaMax1 = std::max( std::abs( beta_dynamicParticle - beta_postStepPoint ), 
                                        std::abs( beta_dynamicParticle - beta_velocity_track ) );
    beta_deltaMax1 =
      std::max( beta_deltaMax1, std::abs( beta_dynamicParticle - beta_velocity_postStepPoint ) );
    const G4double beta_deltaMax2 = std::abs( beta_dynamicParticle - beta_p_over_etot );
    //G4cout << "\t beta_deltaMax{1,2} = " << beta_deltaMax1 << " , " << beta_deltaMax2 << G4endl;
    const G4double beta_val = beta_dynamicParticle;  // To be used later
    // Lorentz gamma of the primary particle at the decay
    const G4double gamma_postStepPoint = theStep->GetPostStepPoint()->GetGamma();
    const G4double gamma_from_e_over_m = etot_val / mass_val;
    const G4double gamma_deltaMax1 = std::abs( gamma_postStepPoint - gamma_from_e_over_m );
    G4double gamma_from_beta = 0.0;
    G4double gamma_deltaMax2 = 0.0;
    G4double gamma_deltaMax3 = 0.0;
    if ( beta_val < 1.0 ) {
      gamma_from_beta = 1.0 / std::sqrt( 1.0 - beta_val*beta_val );
      gamma_deltaMax2 = std::abs( gamma_postStepPoint - gamma_from_beta );
      gamma_deltaMax3 = std::abs( gamma_from_e_over_m - gamma_from_beta );
    }
    const G4double gamma_val = gamma_postStepPoint;  // To be used later;
    //G4cout << "\t gamma_deltaMax{1,2,3} = " << gamma_deltaMax1 << " , " << gamma_deltaMax2
    //       << " , " << gamma_deltaMax3  << " ; gamma_postStepPoint = " << gamma_postStepPoint
    //       << " ; gamma_from_e_over_m = " << gamma_from_e_over_m << " ; gamma_from_beta = "
    //       << gamma_from_beta << G4endl;
    // Proper time of the primary particle at the decay
    const G4double t_proper_track         = theStep->GetTrack()->GetProperTime();
    const G4double t_proper_postStepPoint = theStep->GetPostStepPoint()->GetProperTime();
    const G4double t_proper_deltaMax = std::abs( t_proper_track - t_proper_postStepPoint );
    //G4cout << "\t t_proper_deltaMax [fs] = " << t_proper_deltaMax / femtosecond << G4endl;     
    const G4double t_proper_val = t_proper_track;  // To be used later
    // Lab time of the primary particle at the decay
    // (Note: it would be wrong to trying to compute this lab time from the
    //        above proper time via the simple formula:
    //          const G4double t_lab_from_gamma = t_proper_val * gamma_val;
    //        because the gamma value of the primary particle has changed
    //        during its lifetime.) 
    const G4double t_local_track          = theStep->GetTrack()->GetLocalTime();
    const G4double t_local_postStepPoint  = theStep->GetPostStepPoint()->GetLocalTime();
    const G4double t_global_track         = theStep->GetTrack()->GetGlobalTime();
    const G4double t_global_postStepPoint = theStep->GetPostStepPoint()->GetGlobalTime();
    G4double t_lab_deltaMax = std::max( std::abs( t_local_track - t_local_postStepPoint ),
                                        std::abs( t_local_track - t_global_track ) );
    t_lab_deltaMax = std::max( t_lab_deltaMax,
                               std::abs( t_local_track - t_global_postStepPoint ) );
    //G4cout << "\t t_lab_deltaMax [fs] = " << t_lab_deltaMax / femtosecond << G4endl;
    const G4double t_lab_val = t_local_track;  // To be used later
    // "MC-truth" decay radius of the primary particle at the decay
    // (defined as the one that would happen if there are neither magnetic field effects
    // nor interactions with matter).
    const G4double primaryBeta =
      fPrimaryParticleInitialMomentum / fPrimaryParticleInitialTotalEnergy;
    const G4double mc_truth_rPos1 = t_lab_val * fPrimaryParticleInitialBeta * CLHEP::c_light;
    const G4double mc_truth_rPos2 = t_lab_val * primaryBeta * CLHEP::c_light;
    const G4double mc_truth_rPos_deltaMax = std::abs( mc_truth_rPos1 - mc_truth_rPos2 );
    fMeanMc_truth_rPos_deltaMax += mc_truth_rPos_deltaMax;
    //G4cout << "\t mc_truth_rPos_deltaMax [mum] = "
    //       << mc_truth_rPos_deltaMax / CLHEP::micrometer << G4endl;
    if ( mc_truth_rPos_deltaMax > ToleranceDeltaDecayRadius() ) {
      //G4cout << std::setprecision(6)
      //       << " Large : mc_truth_rPos_deltaMax [mum]="
      //       << mc_truth_rPos_deltaMax / CLHEP::micrometer
      //       << " ; " << mc_truth_rPos1 << " , " << mc_truth_rPos2 << " mm" << G4endl;
      if ( fRunPtr ) fRunPtr->IncrementNumber_mc_truth_rPos_deltaMax_above();
    }
    const G4double mc_truth_rPos_val = mc_truth_rPos1;  // To be used later
    // Keep note of the biggest discrepancies
    fMaxEkin_deltaMax          = std::max( fMaxEkin_deltaMax, ekin_deltaMax );
    fMaxEtot_deltaMax          = std::max( fMaxEtot_deltaMax, etot_deltaMax );
    fMaxP_deltaMax             = std::max( fMaxP_deltaMax, p_deltaMax );
    fMaxPdir_deltaMax          = std::max( fMaxPdir_deltaMax, pdir_deltaMax );
    fMaxMass_deltaMax1         = std::max( fMaxMass_deltaMax1, mass_deltaMax1 );
    fMaxMass_deltaMax2         = std::max( fMaxMass_deltaMax2, mass_deltaMax2 );
    fMaxMass_deltaMax3         = std::max( fMaxMass_deltaMax3, mass_deltaMax3 );
    fMaxBeta_deltaMax1         = std::max( fMaxBeta_deltaMax1, beta_deltaMax1 );
    fMaxBeta_deltaMax2         = std::max( fMaxBeta_deltaMax2, beta_deltaMax2 );
    fMaxGamma_deltaMax1        = std::max( fMaxGamma_deltaMax1, gamma_deltaMax1 );
    fMaxGamma_deltaMax2        = std::max( fMaxGamma_deltaMax2, gamma_deltaMax2 );
    fMaxGamma_deltaMax3        = std::max( fMaxGamma_deltaMax3, gamma_deltaMax3 );
    fMaxT_lab_deltaMax         = std::max( fMaxT_lab_deltaMax, t_lab_deltaMax );
    fMaxT_proper_deltaMax      = std::max( fMaxT_proper_deltaMax, t_proper_deltaMax );
    fMaxMc_truth_rPos_deltaMax = std::max( fMaxMc_truth_rPos_deltaMax, mc_truth_rPos_deltaMax );
    //--- End consistency checks ---

    // Global position
    const G4double xPos = theStep->GetPostStepPoint()->GetPosition().x();
    const G4double yPos = theStep->GetPostStepPoint()->GetPosition().y();
    const G4double zPos = theStep->GetPostStepPoint()->GetPosition().z();
    const G4double rPos = std::sqrt( xPos*xPos + yPos*yPos + zPos*zPos );
    // I have verified that for this case in which only primaries are considered, the
    // "GetGlobalTime()" is the same as "GetLocalTime()" (the one we use).
    // Moreover, this value is also the same as "GetProperTime()"*gamma .
    G4double tPos = theStep->GetPostStepPoint()->GetLocalTime();
    // The "MC-truth" decay radius is defined as the one that would happen if there are
    // neither magnetic field effects nor interactions with matter.
    const G4double mc_truth_rPos = tPos * fPrimaryParticleInitialBeta * CLHEP::c_light;
    const G4double rDeltaPos = mc_truth_rPos - rPos;
    const G4double eKin = theStep->GetPostStepPoint()->GetKineticEnergy();
    const G4double xMom = theStep->GetPostStepPoint()->GetMomentum().x();
    const G4double yMom = theStep->GetPostStepPoint()->GetMomentum().y();
    const G4double zMom = theStep->GetPostStepPoint()->GetMomentum().z();
    // The compute here the angular deflection, in degrees, between the initial direction of
    // the primary particle - which is along the x-axis, and its direction when it decays.
    G4double xDirection = std::min( theStep->GetPostStepPoint()->GetMomentumDirection().x(), 1.0 );
    if ( xDirection < -1.0 ) xDirection = -1.0;
    const G4double deflection_angle_in_degrees = 57.29*std::acos( xDirection );
    const G4double delta_ekin = fPrimaryParticleInitialKineticEnergy - eKin;
    //G4cout << std::setprecision(6)
    //       << " Decay: tPos[ns]=" << tPos << " ; rPos[mm]=" << rPos << " ; deltaR[mum]="
    //       << rDeltaPos /CLHEP::micrometer << " ; deltaEkin[MeV]=" << delta_ekin
    //       << " ; deltaAngle(deg)=" << deflection_angle_in_degrees << G4endl;
    // If the absolute difference between the "MC-truth" decay radius and the real one is above a
    // given threshold, then we notify this special situation in the output, with "LARGE_DELTA_R"
    // for post-processing evaluation. Moreover, in this case, if the "MC-truth" decay radius is
    // smaller than the real one, then we count this unexpected occurrence and we further notify
    // this special situation in the output with "***UNEXPECTED***" for post-processing
    // evaluation.
    if ( std::abs( rDeltaPos ) > ToleranceDeltaDecayRadius() ) {
      //G4cout << "\t LARGE_DELTA_R : mc_truth_rPos[mm]=" << mc_truth_rPos
      //       << " ; rPos[mm]=" << rPos;
      if ( rDeltaPos < 0.0 ) {
        //G4cout << "\t ***UNEXPECTED***";
        if ( fRunPtr ) fRunPtr->IncrementNumberUnexpectedDecays();
      }
      //G4cout << G4endl;
    }
    fMeanDeltaR_primaryDecay += rDeltaPos;
    fMinDeltaR_primaryDecay = std::min( fMinDeltaR_primaryDecay, rDeltaPos );
    fMaxDeltaR_primaryDecay = std::max( fMaxDeltaR_primaryDecay, rDeltaPos );
    fMeanR_primaryDecay += rPos;
    fMinR_primaryDecay = std::min( fMinR_primaryDecay, rPos );
    fMaxR_primaryDecay = std::max( fMaxR_primaryDecay, rPos );
    fMeanX_primaryDecay += xPos;
    fMinX_primaryDecay = std::min( fMinX_primaryDecay, xPos );
    fMaxX_primaryDecay = std::max( fMaxX_primaryDecay, xPos );
    fMeanY_primaryDecay += yPos;
    fMinY_primaryDecay = std::min( fMinY_primaryDecay, yPos );
    fMaxY_primaryDecay = std::max( fMaxY_primaryDecay, yPos );
    fMeanZ_primaryDecay += zPos;
    fMinZ_primaryDecay = std::min( fMinZ_primaryDecay, zPos );
    fMaxZ_primaryDecay = std::max( fMaxZ_primaryDecay, zPos );
    fMeanDeltaAngle_primaryDecay += deflection_angle_in_degrees;
    fMinDeltaAngle_primaryDecay = std::min( fMinDeltaAngle_primaryDecay,
                                            deflection_angle_in_degrees );
    fMaxDeltaAngle_primaryDecay = std::max( fMaxDeltaAngle_primaryDecay,
                                            deflection_angle_in_degrees );      
    fMeanDeltaEkin_primaryDecay += delta_ekin;
    fMinDeltaEkin_primaryDecay = std::min( fMinDeltaEkin_primaryDecay, delta_ekin );
    fMaxDeltaEkin_primaryDecay = std::max( fMaxDeltaEkin_primaryDecay, delta_ekin );
    fMeanEkin_primaryDecay += eKin;
    fMinEkin_primaryDecay = std::min( fMinEkin_primaryDecay, eKin );
    fMaxEkin_primaryDecay = std::max( fMaxEkin_primaryDecay, eKin );
    fMeanPx_primaryDecay += xMom;
    fMinPx_primaryDecay = std::min( fMinPx_primaryDecay, xMom );
    fMaxPx_primaryDecay = std::max( fMaxPx_primaryDecay, xMom );
    fMeanPy_primaryDecay += yMom;
    fMinPy_primaryDecay = std::min( fMinPy_primaryDecay, yMom );
    fMaxPy_primaryDecay = std::max( fMaxPy_primaryDecay, yMom );
    fMeanPz_primaryDecay += zMom;
    fMinPz_primaryDecay = std::min( fMinPz_primaryDecay, zMom );
    fMaxPz_primaryDecay = std::max( fMaxPz_primaryDecay, zMom );

    //--- Extra checks ---
    // Compute the "MC-truth" decay radius using the proper time of the primary particle when
    // it decays.
    // To do this, we would need an "effective" or "average" Lorentz beta and gamma of the primary
    // particle during its lifetime, whereas in practice we have only the Lorentz beta and gamma
    // values at the beginning and at the end when it decays. So, we can get only either an
    // overestimate of the "MC-truth" decay radius - by using the initial Lorentz beta and gamma -
    // or an underestimate of it - by using the Lorentz beta and gamma at the decay.
    // We want to check the average values and the largest values of these wrong estimates.
    const G4double underestimated_mc_truth_rPos =
      t_proper_val * gamma_val * beta_val * CLHEP::c_light;
    const G4double overestimated_mc_truth_rPos =
      t_proper_val * fPrimaryParticleInitialGamma * fPrimaryParticleInitialBeta * CLHEP::c_light;
    const G4double underestimated_mc_truth_rPos_delta =
      underestimated_mc_truth_rPos - mc_truth_rPos_val;
    const G4double overestimated_mc_truth_rPos_delta =
      overestimated_mc_truth_rPos  - mc_truth_rPos_val;
    fMeanUnderestimated_mc_truth_rPos_delta += underestimated_mc_truth_rPos_delta;
    fMeanOverestimated_mc_truth_rPos_delta  += overestimated_mc_truth_rPos_delta;
    //G4cout << "\t underestimated_mc_truth_rPos_delta [mum] = "
    //       << underestimated_mc_truth_rPos_delta / CLHEP::micrometer
    //       << " ; overestimated_mc_truth_rPos_delta [mum] = "
    //       << overestimated_mc_truth_rPos_delta / CLHEP::micrometer << G4endl;
    if ( -underestimated_mc_truth_rPos_delta > ToleranceDeltaDecayRadius() ) {
      //G4cout << std::setprecision(6)
      //       << " Large : underestimated_mc_truth_rPos_delta [mum]="
      //       << underestimated_mc_truth_rPos_delta / CLHEP::micrometer
      //       << " ; " << underestimated_mc_truth_rPos << " , "
      //       << mc_truth_rPos_val << " mm" << G4endl;
      if ( fRunPtr ) fRunPtr->IncrementNumber_underestimated_mc_truth_rPos_delta_above();
    }
    if ( overestimated_mc_truth_rPos_delta > ToleranceDeltaDecayRadius() ) {
      //G4cout << std::setprecision(6)
      //       << " Large : overestimated_mc_truth_rPos_delta [mum]="
      //       << overestimated_mc_truth_rPos_delta / CLHEP::micrometer
      //       << " ; " << overestimated_mc_truth_rPos << " , "
      //       << mc_truth_rPos_val << " mm" << G4endl;
      if ( fRunPtr ) fRunPtr->IncrementNumber_overestimated_mc_truth_rPos_delta_above();
    }    
    const G4double underestimated_rDeltaPos = underestimated_mc_truth_rPos - rPos;
    const G4double overestimated_rDeltaPos  = overestimated_mc_truth_rPos  - rPos;
    fMeanUnderestimated_rDeltaPos += underestimated_rDeltaPos;
    fMeanOverestimated_rDeltaPos  += overestimated_rDeltaPos;
    //G4cout << std::setprecision(6)
    //       << "\t underestimated_rDeltaPos=" << underestimated_rDeltaPos/CLHEP::micrometer
    //       << " ; overestimated_rDeltaPos=" << overestimated_rDeltaPos/CLHEP::micrometer
    //       << " mum" << G4endl;
    if ( -underestimated_rDeltaPos > ToleranceDeltaDecayRadius() ) {
      if ( fRunPtr ) fRunPtr->IncrementNumberLargeUnderestimates();
    }
    if (  overestimated_rDeltaPos  > ToleranceDeltaDecayRadius() ) {
      if ( fRunPtr ) fRunPtr->IncrementNumberLargeOverestimates();
    }
    // Keep note of the biggest discrepancies
    fMinUnderestimated_mc_truth_rPos_delta = std::min( fMinUnderestimated_mc_truth_rPos_delta,
                                                       underestimated_mc_truth_rPos_delta );
    fMaxOverestimated_mc_truth_rPos_delta = std::max( fMaxOverestimated_mc_truth_rPos_delta,
                                                      overestimated_mc_truth_rPos_delta );  
    fMinUnderestimated_rDeltaPos = std::min( fMinUnderestimated_rDeltaPos,
                                             underestimated_rDeltaPos );
    fMaxOverestimated_rDeltaPos = std::max( fMaxOverestimated_rDeltaPos,
                                            overestimated_rDeltaPos );
    //--- End extra checks ---

    // Check numerical errors due to the use of  float  instead of  double : try out
    // several, equivalent computations, taking the one with the largest numerical error.
    const G4float float_xPos =
      static_cast< G4float >( theStep->GetPostStepPoint()->GetPosition().x() );
    const G4float float_yPos =
      static_cast< G4float >( theStep->GetPostStepPoint()->GetPosition().y() );
    const G4float float_zPos =
      static_cast< G4float >( theStep->GetPostStepPoint()->GetPosition().z() );
    const G4float float_rPos =
      std::sqrt( float_xPos*float_xPos + float_yPos*float_yPos + float_zPos*float_zPos );
    const G4float float_tPos =
      static_cast< G4float >( theStep->GetPostStepPoint()->GetLocalTime() );
    const G4float float_initialBeta1 = static_cast< G4float >( fPrimaryParticleInitialBeta );
    const G4float float_initialBeta2 =
      static_cast< G4float >( fPrimaryParticleInitialMomentum ) /
      static_cast< G4float >( fPrimaryParticleInitialTotalEnergy );
    const G4float float_initialGamma = static_cast< G4float >( fPrimaryParticleInitialGamma );
    const G4float float_initialBeta3 =
      std::sqrt( float_initialGamma*float_initialGamma - 1.0 ) / float_initialGamma;
    const G4float float_c_light = static_cast< G4float >( CLHEP::c_light );
    const G4float float_mc_truth_rPos1 = float_tPos * float_initialBeta1 * float_c_light;
    const G4float float_mc_truth_rPos2 = float_tPos * float_initialBeta2 * float_c_light;
    const G4float float_mc_truth_rPos3 = float_tPos * float_initialBeta3 * float_c_light;
    const G4float float_rDeltaPos_0 = static_cast< G4float >( rDeltaPos );
    const G4float float_rDeltaPos_1 = float_mc_truth_rPos1 - float_rPos;
    const G4float float_rDeltaPos_2 = float_mc_truth_rPos2 - float_rPos;
    const G4float float_rDeltaPos_3 = float_mc_truth_rPos3 - float_rPos;
    const G4float float_rDeltaPos_4 = static_cast< G4float >( mc_truth_rPos ) - float_rPos;    
    const G4float float_rDeltaPos_5 = float_mc_truth_rPos1 - static_cast< G4float >( rPos );
    const G4float float_rDeltaPos_6 = float_mc_truth_rPos2 - static_cast< G4float >( rPos );
    const G4float float_rDeltaPos_7 = float_mc_truth_rPos3 - static_cast< G4float >( rPos );
    G4double rDeltaPos_deltaMax = std::max( std::abs( float_rDeltaPos_0 - rDeltaPos ),
                                            std::abs( float_rDeltaPos_1 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_2 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_3 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_4 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_5 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_6 - rDeltaPos ) );
    rDeltaPos_deltaMax = std::max( rDeltaPos_deltaMax,
                                   std::abs( float_rDeltaPos_7 - rDeltaPos ) );
    //G4cout << std::setprecision(6) << " rDeltaPos_deltaMax[mum]=" 
    //       << rDeltaPos_deltaMax / CLHEP::micrometer << G4endl;
    fMaxFloat_rDeltaPos_deltaMax = std::max( fMaxFloat_rDeltaPos_deltaMax, rDeltaPos_deltaMax );
    
    // Get properties of the decay products and check the energy-momentum conservation of the
    // decay
    std::size_t nSec = theStep->GetNumberOfSecondariesInCurrentStep();
    const std::vector< const G4Track* >* ptrVecSecondaries = theStep->GetSecondaryInCurrentStep();
    G4double deltaE = 0.0, deltaPx = 0.0, deltaPy = 0.0, deltaPz = 0.0;
    if ( nSec > 0  &&  ptrVecSecondaries != nullptr ) {
      G4double sumEsecondaries = 0.0;
      G4ThreeVector sumPsecondaries( 0.0, 0.0, 0.0 );
      for ( std::size_t i = 0; i < nSec; ++i ) {
        if ( (*ptrVecSecondaries)[i] ) {
          sumEsecondaries += (*ptrVecSecondaries)[i]->GetTotalEnergy();
          sumPsecondaries += (*ptrVecSecondaries)[i]->GetMomentum();
        }
      }
      deltaE = sumEsecondaries - theStep->GetPostStepPoint()->GetTotalEnergy();
      fMeanViolationE_primaryDecay += deltaE;
      fMinViolationE_primaryDecay = std::min( fMinViolationE_primaryDecay, deltaE );      
      fMaxViolationE_primaryDecay = std::max( fMaxViolationE_primaryDecay, deltaE );      
      if ( std::abs( deltaE ) > ToleranceEPviolations() ) {
        if ( fRunPtr ) fRunPtr->IncrementNumberEviolations();
      }
      deltaPx = sumPsecondaries.x() - xMom;
      fMeanViolationPx_primaryDecay += deltaPx;
      fMinViolationPx_primaryDecay = std::min( fMinViolationPx_primaryDecay, deltaPx );      
      fMaxViolationPx_primaryDecay = std::max( fMaxViolationPx_primaryDecay, deltaPx );
      deltaPy = sumPsecondaries.y() - yMom;
      fMeanViolationPy_primaryDecay += deltaPy;
      fMinViolationPy_primaryDecay = std::min( fMinViolationPy_primaryDecay, deltaPy );
      fMaxViolationPy_primaryDecay = std::max( fMaxViolationPy_primaryDecay, deltaPy );
      deltaPz = sumPsecondaries.z() - zMom;
      fMeanViolationPz_primaryDecay += deltaPz;
      fMinViolationPz_primaryDecay = std::min( fMinViolationPz_primaryDecay, deltaPz );
      fMaxViolationPz_primaryDecay = std::max( fMaxViolationPz_primaryDecay, deltaPz );
      if ( std::abs( deltaPx ) > ToleranceEPviolations()  ||
           std::abs( deltaPy ) > ToleranceEPviolations()  ||
           std::abs( deltaPz ) > ToleranceEPviolations() ) {
        if ( fRunPtr ) fRunPtr->IncrementNumberPviolations();
      }
    } else {
      if ( fRunPtr ) fRunPtr->IncrementNumberBadPrimaryDecays();
    }

    if ( fRunPtr ) {
      fRunPtr->IncrementNumberDecays();
      fRunPtr->SetDecayT( tPos );
      fRunPtr->SetDecayR_mc_truth( mc_truth_rPos );
      fRunPtr->SetDecayR( rPos );
      fRunPtr->SetDecayX( xPos );
      fRunPtr->SetDecayY( yPos );
      fRunPtr->SetDecayZ( zPos );
      fRunPtr->SetDeltaDecayR( rDeltaPos );
      fRunPtr->SetDeflectionAngle( deflection_angle_in_degrees );
      fRunPtr->SetDeltaEkin( delta_ekin );
      fRunPtr->SetDecayEkin( eKin );
      fRunPtr->SetDecayPx( xMom );
      fRunPtr->SetDecayPy( yMom );
      fRunPtr->SetDecayPz( zMom );
      fRunPtr->SetDecayEtotViolation( deltaE );
      fRunPtr->SetDecayPxViolation( deltaPx );
      fRunPtr->SetDecayPyViolation( deltaPy );
      fRunPtr->SetDecayPzViolation( deltaPz );
      fRunPtr->SetMaxEkin_deltaMax( ekin_deltaMax );
      fRunPtr->SetMaxEtot_deltaMax( etot_deltaMax );
      fRunPtr->SetMaxP_deltaMax( p_deltaMax );
      fRunPtr->SetMaxPdir_deltaMax( pdir_deltaMax );
      fRunPtr->SetMaxMass_deltaMax1( mass_deltaMax1 );
      fRunPtr->SetMaxMass_deltaMax2( mass_deltaMax2 );
      fRunPtr->SetMaxMass_deltaMax3( mass_deltaMax3 );
      fRunPtr->SetMaxBeta_deltaMax1( beta_deltaMax1 );
      fRunPtr->SetMaxBeta_deltaMax2( beta_deltaMax2 );
      fRunPtr->SetMaxGamma_deltaMax1( gamma_deltaMax1 );
      fRunPtr->SetMaxGamma_deltaMax2( gamma_deltaMax2 );
      fRunPtr->SetMaxGamma_deltaMax3( gamma_deltaMax3 );
      fRunPtr->SetMaxT_proper_deltaMax( t_proper_deltaMax );
      fRunPtr->SetMaxT_lab_deltaMax( t_lab_deltaMax );
      fRunPtr->SetMaxMc_truth_rPos_deltaMax( mc_truth_rPos_deltaMax );
      fRunPtr->SetMinUnderestimated_mc_truth_rPos_delta( underestimated_mc_truth_rPos_delta );
      fRunPtr->SetMaxOverestimated_mc_truth_rPos_delta( overestimated_mc_truth_rPos_delta ); 
      fRunPtr->SetMinUnderestimated_rDeltaPos( underestimated_rDeltaPos );
      fRunPtr->SetMaxOverestimated_rDeltaPos( overestimated_rDeltaPos );
      fRunPtr->SetMaxFloat_rDeltaPos_deltaMax( fMaxFloat_rDeltaPos_deltaMax ); 
    }

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
