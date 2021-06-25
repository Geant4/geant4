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

#include "G4AdjointAlongStepWeightCorrection.hh"

#include "G4AdjointCSManager.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"

///////////////////////////////////////////////////////
G4AdjointAlongStepWeightCorrection::G4AdjointAlongStepWeightCorrection(
  const G4String& name, G4ProcessType type)
  : G4VContinuousProcess(name, type)
{
  fParticleChange = new G4ParticleChange();
  fCSManager      = G4AdjointCSManager::GetAdjointCSManager();
}

///////////////////////////////////////////////////////
G4AdjointAlongStepWeightCorrection::~G4AdjointAlongStepWeightCorrection()
{
  delete fParticleChange;
}

///////////////////////////////////////////////////////
void G4AdjointAlongStepWeightCorrection::ProcessDescription(
  std::ostream& out) const
{
  out <<
  "Continuous processes act on adjoint particles to continuously correct their "
  "weight during the adjoint reverse tracking. This process is needed when "
  "the adjoint cross sections are not scaled such that the total adjoint cross "
  "section matches the total forward cross section. By default the mode where "
  "the total adjoint cross section is equal to the total forward cross section "
  "is used and therefore this along step weightcorrection factor is 1. However "
  "in some cases (some energy ranges) the total forward cross section or the "
  "total adjoint cross section can be zero. In this case the along step weight "
  "correction is needed and is given by exp(-(Sigma_tot_adj-Sigma_tot_fwd).dx)"
  "\n";
}

///////////////////////////////////////////////////////
G4VParticleChange* G4AdjointAlongStepWeightCorrection::AlongStepDoIt(
  const G4Track& track, const G4Step& step)
{
  fParticleChange->Initialize(track);

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();

  G4double Tkin = step.GetPostStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* thePartDef = const_cast<G4ParticleDefinition*>(
    track.GetDynamicParticle()->GetDefinition());
  G4double weight_correction = fCSManager->GetContinuousWeightCorrection(
    thePartDef, fPreStepKinEnergy, Tkin, fCurrentCouple, length);

  // Caution!!!
  // It is important to select the weight of the post_step_point as the current
  // weight and not the weight of the track, as the weight of the track is
  // changed after having applied all the along_step_do_it.
  G4double new_weight =
    weight_correction * step.GetPostStepPoint()->GetWeight();

  // The following test check for zero weight.
  // This happens after weight correction of gamma for photo electric effect.
  // When the new weight is 0 it will be later on considered as NaN by G4.
  // Therefore we put a lower limit of 1.e-300. for new_weight
  if(new_weight == 0. || (new_weight <= 0. && new_weight > 0.))
  {
    new_weight = 1.e-300;
  }

  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);

  return fParticleChange;
}

///////////////////////////////////////////////////////
G4double G4AdjointAlongStepWeightCorrection::GetContinuousStepLimit(
  const G4Track& track, G4double, G4double, G4double&)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  fPreStepKinEnergy = track.GetKineticEnergy();
  return DBL_MAX;
}
