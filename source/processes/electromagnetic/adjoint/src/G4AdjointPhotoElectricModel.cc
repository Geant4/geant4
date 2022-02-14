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

#include "G4AdjointPhotoElectricModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Gamma.hh"
#include "G4ParticleChange.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4TrackStatus.hh"

////////////////////////////////////////////////////////////////////////////////
G4AdjointPhotoElectricModel::G4AdjointPhotoElectricModel()
  : G4VEmAdjointModel("AdjointPEEffect")

{
  SetUseMatrix(false);
  SetApplyCutInRange(false);

  fAdjEquivDirectPrimPart   = G4AdjointGamma::AdjointGamma();
  fAdjEquivDirectSecondPart = G4AdjointElectron::AdjointElectron();
  fDirectPrimaryPart        = G4Gamma::Gamma();
  fSecondPartSameType       = false;
  fDirectModel              = new G4PEEffectFluoModel();
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointPhotoElectricModel::~G4AdjointPhotoElectricModel() {}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointPhotoElectricModel::SampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  if(isScatProjToProj)
    return;

  // Compute the fTotAdjointCS vectors if not already done for the current
  // couple and electron energy
  const G4DynamicParticle* aDynPart = aTrack.GetDynamicParticle();
  G4double electronEnergy           = aDynPart->GetKineticEnergy();
  G4ThreeVector electronDirection   = aDynPart->GetMomentumDirection();
  fPreStepAdjointCS =
    fTotAdjointCS;  // The last computed CS was at pre step point
  AdjointCrossSection(aTrack.GetMaterialCutsCouple(), electronEnergy,
                      isScatProjToProj);
  fPostStepAdjointCS = fTotAdjointCS;

  // Sample element
  const G4ElementVector* theElementVector =
    fCurrentMaterial->GetElementVector();
  size_t nelm      = fCurrentMaterial->GetNumberOfElements();
  G4double rand_CS = G4UniformRand() * fXsec[nelm - 1];
  for(fIndexElement = 0; fIndexElement < nelm - 1; ++fIndexElement)
  {
    if(rand_CS < fXsec[fIndexElement])
      break;
  }

  // Sample shell and binding energy
  G4int nShells = (*theElementVector)[fIndexElement]->GetNbOfAtomicShells();
  rand_CS       = fShellProb[fIndexElement][nShells - 1] * G4UniformRand();
  G4int i;
  for(i = 0; i < nShells - 1; ++i)
  {
    if(rand_CS < fShellProb[fIndexElement][i])
      break;
  }
  G4double gammaEnergy =
    electronEnergy + (*theElementVector)[fIndexElement]->GetAtomicShell(i);

  // Sample cos theta
  // Copy of the G4PEEfectFluoModel cos theta sampling method
  // ElecCosThetaDistribution. This method cannot be used directly from
  // G4PEEffectFluoModel because it is a friend method.
  G4double cos_theta = 1.;
  G4double gamma     = 1. + electronEnergy / electron_mass_c2;
  if(gamma <= 5.)
  {
    G4double beta = std::sqrt(gamma * gamma - 1.) / gamma;
    G4double b    = 0.5 * gamma * (gamma - 1.) * (gamma - 2.);

    G4double rndm, term, greject, grejsup;
    if(gamma < 2.)
      grejsup = gamma * gamma * (1. + b - beta * b);
    else
      grejsup = gamma * gamma * (1. + b + beta * b);

    do
    {
      rndm      = 1. - 2. * G4UniformRand();
      cos_theta = (rndm + beta) / (rndm * beta + 1.);
      term      = 1. - beta * cos_theta;
      greject = (1. - cos_theta * cos_theta) * (1. + b * term) / (term * term);
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while(greject < G4UniformRand() * grejsup);
  }

  // direction of the adjoint gamma electron
  G4double sin_theta = std::sqrt(1. - cos_theta * cos_theta);
  G4double phi       = twopi * G4UniformRand();
  G4double dirx      = sin_theta * std::cos(phi);
  G4double diry      = sin_theta * std::sin(phi);
  G4double dirz      = cos_theta;
  G4ThreeVector adjoint_gammaDirection(dirx, diry, dirz);
  adjoint_gammaDirection.rotateUz(electronDirection);

  // Weight correction
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), electronEnergy,
                        gammaEnergy, isScatProjToProj);

  // Create secondary and modify fParticleChange
  G4DynamicParticle* anAdjointGamma = new G4DynamicParticle(
    G4AdjointGamma::AdjointGamma(), adjoint_gammaDirection, gammaEnergy);

  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->AddSecondary(anAdjointGamma);
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointPhotoElectricModel::CorrectPostStepWeight(
  G4ParticleChange* fParticleChange, G4double old_weight,
  G4double adjointPrimKinEnergy, G4double projectileKinEnergy, G4bool)
{
  G4double new_weight = old_weight;

  G4double w_corr =
    G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection() /
    fFactorCSBiasing;
  w_corr *= fPostStepAdjointCS / fPreStepAdjointCS;

  new_weight *= w_corr * projectileKinEnergy / adjointPrimKinEnergy;
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointPhotoElectricModel::AdjointCrossSection(
  const G4MaterialCutsCouple* aCouple, G4double electronEnergy,
  G4bool isScatProjToProj)
{
  if(isScatProjToProj)
    return 0.;

  G4double totBiasedAdjointCS = 0.;
  if(aCouple != fCurrentCouple || fCurrenteEnergy != electronEnergy)
  {
    fTotAdjointCS = 0.;
    DefineCurrentMaterialAndElectronEnergy(aCouple, electronEnergy);
    const G4ElementVector* theElementVector =
      fCurrentMaterial->GetElementVector();
    const G4double* theAtomNumDensityVector =
      fCurrentMaterial->GetVecNbOfAtomsPerVolume();
    size_t nelm = fCurrentMaterial->GetNumberOfElements();
    for(fIndexElement = 0; fIndexElement < nelm; ++fIndexElement)
    {
      fTotAdjointCS += AdjointCrossSectionPerAtom(
                         (*theElementVector)[fIndexElement], electronEnergy) *
                       theAtomNumDensityVector[fIndexElement];
      fXsec[fIndexElement] = fTotAdjointCS;
    }

    totBiasedAdjointCS = std::min(fTotAdjointCS, 0.01);
    fFactorCSBiasing   = totBiasedAdjointCS / fTotAdjointCS;
  }
  return totBiasedAdjointCS;
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointPhotoElectricModel::AdjointCrossSectionPerAtom(
  const G4Element* anElement, G4double electronEnergy)
{
  G4int nShells        = anElement->GetNbOfAtomicShells();
  G4double Z           = anElement->GetZ();
  G4double gammaEnergy = electronEnergy + anElement->GetAtomicShell(0);
  G4double CS          = fDirectModel->ComputeCrossSectionPerAtom(
    G4Gamma::Gamma(), gammaEnergy, Z, 0., 0., 0.);
  G4double adjointCS = 0.;
  if(CS > 0.)
    adjointCS += CS / gammaEnergy;
  fShellProb[fIndexElement][0] = adjointCS;
  for(G4int i = 1; i < nShells; ++i)
  {
    G4double Bi1 = anElement->GetAtomicShell(i - 1);
    G4double Bi  = anElement->GetAtomicShell(i);
    if(electronEnergy < Bi1 - Bi)
    {
      gammaEnergy = electronEnergy + Bi;
      CS          = fDirectModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),
                                                    gammaEnergy, Z, 0., 0., 0.);
      if(CS > 0.)
        adjointCS += CS / gammaEnergy;
    }
    fShellProb[fIndexElement][i] = adjointCS;
  }
  adjointCS *= electronEnergy;
  return adjointCS;
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointPhotoElectricModel::DefineCurrentMaterialAndElectronEnergy(
  const G4MaterialCutsCouple* couple, G4double anEnergy)
{
  fCurrentCouple   = const_cast<G4MaterialCutsCouple*>(couple);
  fCurrentMaterial = const_cast<G4Material*>(couple->GetMaterial());
  fCurrenteEnergy = anEnergy;
  fDirectModel->SetCurrentCouple(couple);
}
