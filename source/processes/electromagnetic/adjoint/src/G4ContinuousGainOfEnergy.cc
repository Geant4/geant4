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

#include "G4ContinuousGainOfEnergy.hh"

#include "G4EmCorrections.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4VEmModel.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VParticleChange.hh"

///////////////////////////////////////////////////////
G4ContinuousGainOfEnergy::G4ContinuousGainOfEnergy(const G4String& name,
                                                   G4ProcessType type)
  : G4VContinuousProcess(name, type)
{}

///////////////////////////////////////////////////////
G4ContinuousGainOfEnergy::~G4ContinuousGainOfEnergy() {}

///////////////////////////////////////////////////////
void G4ContinuousGainOfEnergy::ProcessDescription(std::ostream& out) const
{
  out << "Continuous process acting on adjoint particles to compute the "
         "continuous gain of energy of charged particles when they are "
         "tracked back.\n";
}

///////////////////////////////////////////////////////
void G4ContinuousGainOfEnergy::SetDirectParticle(G4ParticleDefinition* p)
{
  fDirectPartDef = p;
  if(fDirectPartDef->GetParticleType() == "nucleus")
  {
    fIsIon     = true;
    fMassRatio = proton_mass_c2 / fDirectPartDef->GetPDGMass();
  }
}

///////////////////////////////////////////////////////
G4VParticleChange* G4ContinuousGainOfEnergy::AlongStepDoIt(const G4Track& track,
                                                           const G4Step& step)
{
  // Caution in this method the step length should be the true step length
  // A problem is that this is computed by the multiple scattering that does
  // not know the energy at the end of the adjoint step. This energy is used
  // during the forward sim. Nothing we can really do against that at this
  // time. This is inherent to the MS method

  aParticleChange.Initialize(track);

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  G4double degain = 0.0;

  // Compute this for weight change after continuous energy loss
  G4double DEDX_before =
    fDirectEnergyLossProcess->GetDEDX(fPreStepKinEnergy, fCurrentCouple);

  // For the fluctuation we generate a new dynamic particle with energy
  // = preEnergy+egain and then compute the fluctuation given in the direct
  // case.
  G4DynamicParticle* dynParticle = new G4DynamicParticle();
  *dynParticle                   = *(track.GetDynamicParticle());
  dynParticle->SetDefinition(fDirectPartDef);
  G4double Tkin = dynParticle->GetKineticEnergy();

  G4double dlength = length;
  if(Tkin != fPreStepKinEnergy && fIsIon)
  {
    G4double chargeSqRatio = fCurrentModel->GetChargeSquareRatio(
      fDirectPartDef, fCurrentMaterial, Tkin);
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio, chargeSqRatio);
  }

  G4double r = fDirectEnergyLossProcess->GetRange(Tkin, fCurrentCouple);
  if(dlength <= fLinLossLimit * r)
  {
    degain = DEDX_before * dlength;
  }
  else
  {
    G4double x = r + dlength;
    G4double E = fDirectEnergyLossProcess->GetKineticEnergy(x, fCurrentCouple);
    if(fIsIon)
    {
      G4double chargeSqRatio = fCurrentModel->GetChargeSquareRatio(
        fDirectPartDef, fCurrentMaterial, E);
      fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio, chargeSqRatio);
      G4double x1 = fDirectEnergyLossProcess->GetRange(E, fCurrentCouple);

      G4int ii              = 0;
      constexpr G4int iimax = 100;
      while(std::abs(x - x1) > 0.01 * x)
      {
        E = fDirectEnergyLossProcess->GetKineticEnergy(x, fCurrentCouple);
        chargeSqRatio = fCurrentModel->GetChargeSquareRatio(
          fDirectPartDef, fCurrentMaterial, E);
        fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio,
                                                       chargeSqRatio);
        x1 = fDirectEnergyLossProcess->GetRange(E, fCurrentCouple);
        ++ii;
        if(ii >= iimax)
        {
          break;
        }
      }
    }

    degain = E - Tkin;
  }
  G4double tmax = fCurrentModel->MaxSecondaryKinEnergy(dynParticle);
  fCurrentTcut = std::min(fCurrentTcut, tmax);

  dynParticle->SetKineticEnergy(Tkin + degain);

  // Corrections, which cannot be tabulated for ions
  fCurrentModel->CorrectionsAlongStep(fCurrentCouple, dynParticle, dlength, degain);

  // Sample fluctuations
  G4double deltaE = 0.;
  if(fLossFluctuationFlag)
  {
    deltaE = fCurrentModel->GetModelOfFluctuations()->SampleFluctuations(
      fCurrentCouple, dynParticle, fCurrentTcut, tmax, dlength, degain) 
      - degain;
  }

  G4double egain = degain + deltaE;
  if(egain <= 0.)
    egain = degain;
  Tkin += egain;
  dynParticle->SetKineticEnergy(Tkin);

  delete dynParticle;

  if(fIsIon)
  {
    G4double chargeSqRatio = fCurrentModel->GetChargeSquareRatio(
      fDirectPartDef, fCurrentMaterial, Tkin);
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio, chargeSqRatio);
  }

  G4double DEDX_after = fDirectEnergyLossProcess->GetDEDX(Tkin, fCurrentCouple);
  G4double weight_correction = DEDX_after / DEDX_before;

  aParticleChange.ProposeEnergy(Tkin);

  // Caution!!! It is important to select the weight of the post_step_point
  // as the current weight and not the weight of the track, as the  weight of
  // the track is changed after having applied all the along_step_do_it.

  G4double new_weight =
    weight_correction * step.GetPostStepPoint()->GetWeight();
  aParticleChange.SetParentWeightByProcess(false);
  aParticleChange.ProposeParentWeight(new_weight);

  return &aParticleChange;
}

///////////////////////////////////////////////////////
void G4ContinuousGainOfEnergy::SetLossFluctuations(G4bool val)
{
  if(val && !fLossFluctuationArePossible)
    return;
  fLossFluctuationFlag = val;
}

///////////////////////////////////////////////////////
G4double G4ContinuousGainOfEnergy::GetContinuousStepLimit(const G4Track& track,
                                                          G4double, G4double,
                                                          G4double&)
{
  DefineMaterial(track.GetMaterialCutsCouple());

  fPreStepKinEnergy = track.GetKineticEnergy();
  fCurrentModel     = fDirectEnergyLossProcess->SelectModelForMaterial(
    track.GetKineticEnergy() * fMassRatio, fCurrentCoupleIndex);
  G4double emax_model           = fCurrentModel->HighEnergyLimit();
  G4double preStepChargeSqRatio = 0.;
  if(fIsIon)
  {
    G4double chargeSqRatio = fCurrentModel->GetChargeSquareRatio(
      fDirectPartDef, fCurrentMaterial, fPreStepKinEnergy);
    preStepChargeSqRatio = chargeSqRatio;
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio,
                                                   preStepChargeSqRatio);
  }

  G4double maxE = 1.1 * fPreStepKinEnergy;

  if(fPreStepKinEnergy < fCurrentTcut)
    maxE = std::min(fCurrentTcut, maxE);

  maxE = std::min(emax_model * 1.001, maxE);

  G4double preStepRange =
    fDirectEnergyLossProcess->GetRange(fPreStepKinEnergy, fCurrentCouple);

  if(fIsIon)
  {
    G4double chargeSqRatioAtEmax = fCurrentModel->GetChargeSquareRatio(
      fDirectPartDef, fCurrentMaterial, maxE);
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio,
                                                   chargeSqRatioAtEmax);
  }

  G4double r1 = fDirectEnergyLossProcess->GetRange(maxE, fCurrentCouple);

  if(fIsIon)
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio,
                                                   preStepChargeSqRatio);

  return std::max(r1 - preStepRange, 0.001 * mm);
}

///////////////////////////////////////////////////////
void G4ContinuousGainOfEnergy::SetDynamicMassCharge(const G4Track&,
                                                    G4double energy)
{
  G4double ChargeSqRatio =
    G4LossTableManager::Instance()->EmCorrections()->EffectiveChargeSquareRatio(
      fDirectPartDef, fCurrentMaterial, energy);
  if(fDirectEnergyLossProcess)
    fDirectEnergyLossProcess->SetDynamicMassCharge(fMassRatio, ChargeSqRatio);
}
