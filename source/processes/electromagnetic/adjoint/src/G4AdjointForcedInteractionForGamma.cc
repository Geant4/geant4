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

#include "G4AdjointForcedInteractionForGamma.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointGamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmAdjointModel.hh"

G4AdjointForcedInteractionForGamma::G4AdjointForcedInteractionForGamma(
  G4String process_name)
  : G4VContinuousDiscreteProcess(process_name)
  , fAdjointComptonModel(nullptr)
  , fAdjointBremModel(nullptr)
{
  fCSManager      = G4AdjointCSManager::GetAdjointCSManager();
  fParticleChange = new G4ParticleChange();
}

//////////////////////////////////////////////////////////////////////////////
G4AdjointForcedInteractionForGamma::~G4AdjointForcedInteractionForGamma()
{
  if(fParticleChange)
    delete fParticleChange;
}

//////////////////////////////////////////////////////////////////////////////
void G4AdjointForcedInteractionForGamma::ProcessDescription(
  std::ostream& out) const
{
  out << "Forced interaction for gamma.\n";
}

//////////////////////////////////////////////////////////////////////////////
void G4AdjointForcedInteractionForGamma::BuildPhysicsTable(
  const G4ParticleDefinition&)
{
  fCSManager->BuildCrossSectionMatrices();  // it will be done just once
  fCSManager->BuildTotalSigmaTables();
}

// Note on weight correction for forced interaction.
// For the forced interaction applied here we use a truncated exponential law
// for the probability of survival over a fixed total length. This is done by
// using a linear transformation of the non-biased probability survival. In
// math this is written P'(x)=C1P(x)+C2 , with P(x)=exp(-sum(sigma_ixi)) . x and
// L can cross different volumes with different cross section sigma. For forced
// interaction, we get the limit conditions:
// P'(L)=0  and  P'(0)=1       (L can be used over different volumes)
// From simple solving of linear equations we
// get C1=1/(1-P(L)) and C2=-P(L)/(1-P(L))
// P'(x)=(P(x)-P(L))/(1-P(L))
// For the probability over a step x1 to x2, P'(x1->x2)=P'(x2)/P'(x1).
// The effective cross
// section is defined -d(P'(x))/dx/P'(x).
// We get therefore
// sigma_eff = C1sigmaP(x)/(C1P(x)+C2) = sigmaP(x)/(P(x)+C2/C1)
//           = sigmaP(x)/(P(x)-P(L)) = sigma/(1-P(L)/P(x))
//////////////////////////////////////////////////////////////////////////////
G4VParticleChange* G4AdjointForcedInteractionForGamma::PostStepDoIt(
  const G4Track& track, const G4Step&)
{
  fParticleChange->Initialize(track);
  // For the free flight gamma no interaction occurs but a gamma with same
  // properties is produced for further forced interaction. It is done at the
  // very beginning of the track so that the weight can be the same
  if(fCopyGammaForForced)
  {
    G4ThreeVector theGammaMomentum = track.GetMomentum();
    fParticleChange->AddSecondary(
      new G4DynamicParticle(G4AdjointGamma::AdjointGamma(), theGammaMomentum));
    fParticleChange->SetParentWeightByProcess(false);
    fParticleChange->SetSecondaryWeightByProcess(false);
  }
  else
  {  // Occurrence of forced interaction
    // Selection of the model to be called
    G4VEmAdjointModel* theSelectedModel = nullptr;
    G4bool is_scat_proj_to_proj_case    = false;
    G4double factor=1.;
    if(!fAdjointComptonModel && !fAdjointBremModel)
      return fParticleChange;
    if(!fAdjointComptonModel)
    {
      theSelectedModel          = fAdjointBremModel;
      is_scat_proj_to_proj_case = false;
      // This is needed because the results of it will be used in the post step
      // do it weight correction inside the model
      fAdjointBremModel->AdjointCrossSection(track.GetMaterialCutsCouple(),
                                             track.GetKineticEnergy(), false);
    }
    else if(!fAdjointBremModel)
    {
      theSelectedModel          = fAdjointComptonModel;
      is_scat_proj_to_proj_case = true;
    }
    else
    {  // Choose the model according to a 50-50 % probability
      G4double bremAdjCS = fAdjointBremModel->AdjointCrossSection(
        track.GetMaterialCutsCouple(), track.GetKineticEnergy(), false);
      if(G4UniformRand()  < 0.5)
      {
        theSelectedModel          = fAdjointBremModel;
        is_scat_proj_to_proj_case = false;
        factor=bremAdjCS/fLastAdjCS/0.5;
      }
      else
      {
        theSelectedModel          = fAdjointComptonModel;
        is_scat_proj_to_proj_case = true;
        factor=(fLastAdjCS-bremAdjCS)/fLastAdjCS/0.5;
      }
    }

    // Compute the weight correction factor
    G4double invEffectiveAdjointCS =
      (1. - std::exp(fNbAdjIntLength - fTotNbAdjIntLength)) / fLastAdjCS/fCSBias;

    // Call the  selected model without correction of the weight in the model
    theSelectedModel->SetCorrectWeightForPostStepInModel(false);
    theSelectedModel
      ->SetAdditionalWeightCorrectionFactorForPostStepOutsideModel(
        factor*fLastAdjCS * invEffectiveAdjointCS);
    theSelectedModel->SampleSecondaries(track, is_scat_proj_to_proj_case,
                                        fParticleChange);
    theSelectedModel->SetCorrectWeightForPostStepInModel(true);

    fContinueGammaAsNewFreeFlight = true;
  }
  return fParticleChange;
}

//////////////////////////////////////////////////////////////////////////////
G4VParticleChange* G4AdjointForcedInteractionForGamma::AlongStepDoIt(
  const G4Track& track, const G4Step&)
{
  fParticleChange->Initialize(track);
  // Compute nb of interactions length over step length
  G4ThreeVector position = track.GetPosition();
  G4double stepLength    = track.GetStep()->GetStepLength();
  G4double ekin          = track.GetKineticEnergy();
  fLastAdjCS = fCSManager->GetTotalAdjointCS(track.GetDefinition(), ekin,
                                             track.GetMaterialCutsCouple());
  G4double nb_fwd_interaction_length_over_step =
    stepLength * fCSManager->GetTotalForwardCS(G4AdjointGamma::AdjointGamma(),
                                               ekin,
                                               track.GetMaterialCutsCouple());

  G4double nb_adj_interaction_length_over_step = stepLength * fLastAdjCS;
  G4double fwd_survival_probability =
    std::exp(-nb_fwd_interaction_length_over_step);
  G4double mc_induced_survival_probability = 1.;

  if(fFreeFlightGamma)
  {  // for free_flight survival probability stays 1
    // Accumulate the number of interaction lengths during free flight of gamma
    fTotNbAdjIntLength += nb_adj_interaction_length_over_step;
    fAccTrackLength += stepLength;
  }
  else
  {
    G4double previous_acc_nb_adj_interaction_length = fNbAdjIntLength;
    fNbAdjIntLength += fCSBias*nb_adj_interaction_length_over_step;
    theNumberOfInteractionLengthLeft -= fCSBias*nb_adj_interaction_length_over_step;

    // protection against rare race condition
    if(std::abs(fTotNbAdjIntLength - previous_acc_nb_adj_interaction_length) <=
       1.e-15)
    {
      mc_induced_survival_probability = 1.e50;
    }
    else
    {
      mc_induced_survival_probability =
        std::exp(-fNbAdjIntLength) - std::exp(-fTotNbAdjIntLength);
      mc_induced_survival_probability /=
        (std::exp(-previous_acc_nb_adj_interaction_length) -
         std::exp(-fTotNbAdjIntLength));
    }
  }
  G4double weight_correction =
    fwd_survival_probability / mc_induced_survival_probability;

  // Caution!!!
  // It is important to select the weight of the post_step_point as the
  // current weight and not the weight of the track, as the weight of the track
  // is changed after having applied all the along_step_do_it.
  G4double new_weight =
    weight_correction * track.GetStep()->GetPostStepPoint()->GetWeight();

  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);

  return fParticleChange;
}

//////////////////////////////////////////////////////////////////////////////
G4double
G4AdjointForcedInteractionForGamma::PostStepGetPhysicalInteractionLength(
  const G4Track& track, G4double, G4ForceCondition* condition)
{
  static G4int lastFreeFlightTrackId = 1000;
  G4int step_id                      = track.GetCurrentStepNumber();
  *condition                         = NotForced;
  fCopyGammaForForced                = false;
  G4int track_id                     = track.GetTrackID();
  fFreeFlightGamma =
    (track_id != lastFreeFlightTrackId + 1 || fContinueGammaAsNewFreeFlight);
  if(fFreeFlightGamma)
  {
    if(step_id == 1 || fContinueGammaAsNewFreeFlight)
    {
      *condition = Forced;
      // A gamma with same conditions will be generate at next post_step do it
      // for the forced interaction
      fCopyGammaForForced           = true;
      lastFreeFlightTrackId         = track_id;
      fAccTrackLength               = 0.;
      fTotNbAdjIntLength            = 0.;
      fContinueGammaAsNewFreeFlight = false;
      return 1.e-90;
    }
    else
    {
      return DBL_MAX;
    }
  }
  else
  {  // compute the interaction length for forced interaction
    if(step_id == 1)
    {
      fCSBias=0.000001/fTotNbAdjIntLength;
      fTotNbAdjIntLength*=fCSBias;
      G4double min_val = std::exp(-fTotNbAdjIntLength);
      theNumberOfInteractionLengthLeft =
        -std::log(min_val + G4UniformRand() * (1. - min_val));
      theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft;
      fNbAdjIntLength                     = 0.;
    }
    G4VPhysicalVolume* thePostPhysVolume =
      track.GetStep()->GetPreStepPoint()->GetPhysicalVolume();
    G4double ekin   = track.GetKineticEnergy();
    G4double postCS = 0.;
    if(thePostPhysVolume)
    {
      postCS = fCSManager->GetTotalAdjointCS(
        G4AdjointGamma::AdjointGamma(), ekin,
        thePostPhysVolume->GetLogicalVolume()->GetMaterialCutsCouple());
    }
    if(postCS > 0.)
      return theNumberOfInteractionLengthLeft / postCS /fCSBias;
    else
      return DBL_MAX;
  }
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointForcedInteractionForGamma::GetContinuousStepLimit(
  const G4Track&, G4double, G4double, G4double&)
{
  return DBL_MAX;
}

////////////////////////////////////////////////////////////////////////////////
// Not used in this process but should be implemented as virtual method
G4double G4AdjointForcedInteractionForGamma::GetMeanFreePath(const G4Track&,
                                                             G4double,
                                                             G4ForceCondition*)
{
  return 0.;
}
