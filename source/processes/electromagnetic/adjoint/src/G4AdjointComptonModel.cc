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

#include "G4AdjointComptonModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Gamma.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"
#include "G4VEmProcess.hh"

////////////////////////////////////////////////////////////////////////////////
G4AdjointComptonModel::G4AdjointComptonModel()
  : G4VEmAdjointModel("AdjointCompton")
{
  SetApplyCutInRange(false);
  SetUseMatrix(false);
  SetUseMatrixPerElement(true);
  SetUseOnlyOneMatrixForAllElements(true);
  fAdjEquivDirectPrimPart   = G4AdjointGamma::AdjointGamma();
  fAdjEquivDirectSecondPart = G4AdjointElectron::AdjointElectron();
  fDirectPrimaryPart        = G4Gamma::Gamma();
  fSecondPartSameType       = false;
  fDirectModel =
    new G4KleinNishinaCompton(G4Gamma::Gamma(), "ComptonDirectModel");
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointComptonModel::~G4AdjointComptonModel() {}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointComptonModel::SampleSecondaries(const G4Track& aTrack,
                                              G4bool isScatProjToProj,
                                              G4ParticleChange* fParticleChange)
{
  if(!fUseMatrix)
    return RapidSampleSecondaries(aTrack, isScatProjToProj, fParticleChange);

  // A recall of the compton scattering law:
  // Egamma2=Egamma1/(1+(Egamma1/E0_electron)(1.-cos_th))
  // Therefore Egamma2_max= Egamma2(cos_th=1) = Egamma1
  // and       Egamma2_min= Egamma2(cos_th=-1) =
  //             Egamma1/(1+2.(Egamma1/E0_electron))
  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();
  G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  // Sample secondary energy
  G4double gammaE1;
  gammaE1 =
    SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy, isScatProjToProj);

  // gammaE2
  G4double gammaE2 = adjointPrimKinEnergy;
  if(!isScatProjToProj)
    gammaE2 = gammaE1 - adjointPrimKinEnergy;

  // Cos th
  G4double cos_th = 1. + electron_mass_c2 * (1. / gammaE1 - 1. / gammaE2);
  if(!isScatProjToProj)
  {
    cos_th =
      (gammaE1 - gammaE2 * cos_th) / theAdjointPrimary->GetTotalMomentum();
  }
  G4double sin_th = 0.;
  if(std::abs(cos_th) > 1.)
  {
    if(cos_th > 0.)
    {
      cos_th = 1.;
    }
    else
      cos_th = -1.;
    sin_th = 0.;
  }
  else
    sin_th = std::sqrt(1. - cos_th * cos_th);

  // gamma0 momentum
  G4ThreeVector dir_parallel = theAdjointPrimary->GetMomentumDirection();
  G4double phi               = G4UniformRand() * twopi;
  G4ThreeVector gammaMomentum1 =
    gammaE1 *
    G4ThreeVector(std::cos(phi) * sin_th, std::sin(phi) * sin_th, cos_th);
  gammaMomentum1.rotateUz(dir_parallel);

  // correct the weight of particles before adding the secondary
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(),
                        adjointPrimKinEnergy, gammaE1, isScatProjToProj);

  if(!isScatProjToProj)
  {  // kill the primary and add a secondary
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->AddSecondary(
      new G4DynamicParticle(fAdjEquivDirectPrimPart, gammaMomentum1));
  }
  else
  {
    fParticleChange->ProposeEnergy(gammaE1);
    fParticleChange->ProposeMomentumDirection(gammaMomentum1.unit());
  }
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointComptonModel::RapidSampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();
  DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());

  G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();

  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  G4double diffCSUsed =
    0.1 * fCurrentMaterial->GetElectronDensity() * twopi_mc2_rcl2;
  G4double gammaE1 = 0.;
  G4double gammaE2 = 0.;
  if(!isScatProjToProj)
  {
    G4double Emax = GetSecondAdjEnergyMaxForProdToProj(adjointPrimKinEnergy);
    G4double Emin = GetSecondAdjEnergyMinForProdToProj(adjointPrimKinEnergy);
    if(Emin >= Emax)
      return;
    G4double f1 = (Emin - adjointPrimKinEnergy) / Emin;
    G4double f2 = (Emax - adjointPrimKinEnergy) / Emax / f1;
    gammaE1 = adjointPrimKinEnergy / (1. - f1 * std::pow(f2, G4UniformRand()));
    gammaE2 = gammaE1 - adjointPrimKinEnergy;
    diffCSUsed =
      diffCSUsed *
      (1. + 2. * std::log(1. + electron_mass_c2 / adjointPrimKinEnergy)) *
      adjointPrimKinEnergy / gammaE1 / gammaE2;
  }
  else
  {
    G4double Emax =
      GetSecondAdjEnergyMaxForScatProjToProj(adjointPrimKinEnergy);
    G4double Emin =
      GetSecondAdjEnergyMinForScatProjToProj(adjointPrimKinEnergy, fTcutSecond);
    if(Emin >= Emax)
      return;
    gammaE2    = adjointPrimKinEnergy;
    gammaE1    = Emin * std::pow(Emax / Emin, G4UniformRand());
    diffCSUsed = diffCSUsed / gammaE1;
  }

  // Weight correction
  // First w_corr is set to the ratio between adjoint total CS and fwd total CS
  G4double w_corr = fOutsideWeightFactor;
  if(fInModelWeightCorr)
  {
    w_corr =
      G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
  }
  // Then another correction is needed because a biased differential CS has
  // been used rather than the one consistent with the direct model

  G4double diffCS =
    DiffCrossSectionPerAtomPrimToScatPrim(gammaE1, gammaE2, 1, 0.);
  if(diffCS > 0.)
    diffCS /= fDirectCS;  // here we have the normalised diffCS
  // And we remultiply by the lambda of the forward process
  diffCS *= fDirectProcess->GetCrossSection(gammaE1, fCurrentCouple);

  w_corr *= diffCS / diffCSUsed;

  G4double new_weight = aTrack.GetWeight() * w_corr;
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);

  G4double cos_th = 1. + electron_mass_c2 * (1. / gammaE1 - 1. / gammaE2);
  if(!isScatProjToProj)
  {
    G4double p_elec = theAdjointPrimary->GetTotalMomentum();
    cos_th          = (gammaE1 - gammaE2 * cos_th) / p_elec;
  }
  G4double sin_th = 0.;
  if(std::abs(cos_th) > 1.)
  {
    if(cos_th > 0.)
    {
      cos_th = 1.;
    }
    else
      cos_th = -1.;
  }
  else
    sin_th = std::sqrt(1. - cos_th * cos_th);

  // gamma0 momentum
  G4ThreeVector dir_parallel = theAdjointPrimary->GetMomentumDirection();
  G4double phi               = G4UniformRand() * twopi;
  G4ThreeVector gammaMomentum1 =
    gammaE1 *
    G4ThreeVector(std::cos(phi) * sin_th, std::sin(phi) * sin_th, cos_th);
  gammaMomentum1.rotateUz(dir_parallel);

  if(!isScatProjToProj)
  {  // kill the primary and add a secondary
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->AddSecondary(
      new G4DynamicParticle(fAdjEquivDirectPrimPart, gammaMomentum1));
  }
  else
  {
    fParticleChange->ProposeEnergy(gammaE1);
    fParticleChange->ProposeMomentumDirection(gammaMomentum1.unit());
  }
}

////////////////////////////////////////////////////////////////////////////////
// The implementation here is correct for energy loss process, for the
// photoelectric and compton scattering the method should be redefined
G4double G4AdjointComptonModel::DiffCrossSectionPerAtomPrimToSecond(
  G4double gamEnergy0, G4double kinEnergyElec, G4double Z, G4double A)
{
  G4double gamEnergy1   = gamEnergy0 - kinEnergyElec;
  G4double dSigmadEprod = 0.;
  if(gamEnergy1 > 0.)
    dSigmadEprod =
      DiffCrossSectionPerAtomPrimToScatPrim(gamEnergy0, gamEnergy1, Z, A);
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointComptonModel::DiffCrossSectionPerAtomPrimToScatPrim(
  G4double gamEnergy0, G4double gamEnergy1, G4double Z, G4double)
{
  // Based on Klein Nishina formula
  // In the forward case (see G4KleinNishinaCompton) the cross section is
  // parametrised but the secondaries are sampled from the Klein Nishina
  // differential cross section. The differential cross section used here
  // is therefore the cross section multiplied by the normalised
  // differential Klein Nishina cross section

  // Klein Nishina Cross Section
  G4double epsilon           = gamEnergy0 / electron_mass_c2;
  G4double one_plus_two_epsi = 1. + 2. * epsilon;
  if(gamEnergy1 > gamEnergy0 || gamEnergy1 < gamEnergy0 / one_plus_two_epsi)
  {
    return 0.;
  }

  G4double CS = std::log(one_plus_two_epsi) *
                (1. - 2. * (1. + epsilon) / (epsilon * epsilon));
  CS +=
    4. / epsilon + 0.5 * (1. - 1. / (one_plus_two_epsi * one_plus_two_epsi));
  CS /= epsilon;
  // Note that the pi*re2*Z factor is neglected because it is suppressed when
  // computing dCS_dE1/CS in the differential cross section

  // Klein Nishina Differential Cross Section
  G4double epsilon1 = gamEnergy1 / electron_mass_c2;
  G4double v        = epsilon1 / epsilon;
  G4double term1    = 1. + 1. / epsilon - 1. / epsilon1;
  G4double dCS_dE1  = 1. / v + v + term1 * term1 - 1.;
  dCS_dE1 *= 1. / epsilon / gamEnergy0;

  // Normalised to the CS used in G4
  fDirectCS = fDirectModel->ComputeCrossSectionPerAtom(
    G4Gamma::Gamma(), gamEnergy0, Z, 0., 0., 0.);

  dCS_dE1 *= fDirectCS / CS;

  return dCS_dE1;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointComptonModel::GetSecondAdjEnergyMaxForScatProjToProj(
  G4double primAdjEnergy)
{
  G4double inv_e_max = 1. / primAdjEnergy - 2. / electron_mass_c2;
  G4double e_max     = GetHighEnergyLimit();
  if(inv_e_max > 0.)
    e_max = std::min(1. / inv_e_max, e_max);
  return e_max;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointComptonModel::GetSecondAdjEnergyMinForProdToProj(
  G4double primAdjEnergy)
{
  G4double half_e = primAdjEnergy / 2.;
  return half_e + std::sqrt(half_e * (electron_mass_c2 + half_e));
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointComptonModel::AdjointCrossSection(
  const G4MaterialCutsCouple* aCouple, G4double primEnergy,
  G4bool isScatProjToProj)
{
  if(fUseMatrix)
    return G4VEmAdjointModel::AdjointCrossSection(aCouple, primEnergy,
                                                  isScatProjToProj);
  DefineCurrentMaterial(aCouple);

  G4float Cross     = 0.;
  G4float Emax_proj = 0.;
  G4float Emin_proj = 0.;
  if(!isScatProjToProj)
  {
    Emax_proj = GetSecondAdjEnergyMaxForProdToProj(primEnergy);
    Emin_proj = GetSecondAdjEnergyMinForProdToProj(primEnergy);
    if(Emax_proj > Emin_proj)
    {
      Cross = 0.1 *
              std::log((Emax_proj - G4float(primEnergy)) * Emin_proj /
                       Emax_proj / (Emin_proj - primEnergy)) *
              (1. + 2. * std::log(G4float(1. + electron_mass_c2 / primEnergy)));
    }
  }
  else
  {
    Emax_proj = GetSecondAdjEnergyMaxForScatProjToProj(primEnergy);
    Emin_proj = GetSecondAdjEnergyMinForScatProjToProj(primEnergy, 0.);
    if(Emax_proj > Emin_proj)
    {
      Cross = 0.1 * std::log(Emax_proj / Emin_proj);
    }
  }

  Cross *= fCurrentMaterial->GetElectronDensity() * twopi_mc2_rcl2;
  fLastCS = Cross;
  return double(Cross);
}
