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

#include "G4AdjointIonIonisationModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggIonModel.hh"
#include "G4GenericIon.hh"
#include "G4NistManager.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"

////////////////////////////////////////////////////////////////////////////////
G4AdjointIonIonisationModel::G4AdjointIonIonisationModel()
  : G4VEmAdjointModel("Adjoint_IonIonisation")
{
  fUseMatrix               = true;
  fUseMatrixPerElement     = true;
  fApplyCutInRange         = true;
  fOneMatrixForAllElements = true;
  fSecondPartSameType      = false;

  // The direct EM Model is taken as BetheBloch. It is only used for the
  // computation of the differential cross section.
  // The Bragg model could be used as an alternative as it offers the same
  // differential cross section

  fBetheBlochDirectEMModel  = new G4BetheBlochModel(G4GenericIon::GenericIon());
  fBraggIonDirectEMModel    = new G4BraggIonModel(G4GenericIon::GenericIon());
  fAdjEquivDirectSecondPart = G4AdjointElectron::AdjointElectron();
  fDirectPrimaryPart        = nullptr;
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointIonIonisationModel::~G4AdjointIonIonisationModel() {}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointIonIonisationModel::SampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();

  // Elastic inverse scattering
  G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
  G4double adjointPrimP         = theAdjointPrimary->GetTotalMomentum();

  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  // Sample secondary energy
  G4double projectileKinEnergy =
    SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy, isScatProjToProj);
  // Caution !!!this weight correction should be always applied
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(),
                        adjointPrimKinEnergy, projectileKinEnergy,
                        isScatProjToProj);

  // Kinematics:
  // we consider a two body elastic scattering for the forward processes where
  // the projectile knock on an e- at rest and gives it part of its  energy
  G4double projectileM0          = fAdjEquivDirectPrimPart->GetPDGMass();
  G4double projectileTotalEnergy = projectileM0 + projectileKinEnergy;
  G4double projectileP2 =
    projectileTotalEnergy * projectileTotalEnergy - projectileM0 * projectileM0;

  // Companion
  G4double companionM0 = fAdjEquivDirectPrimPart->GetPDGMass();
  if(isScatProjToProj)
  {
    companionM0 = fAdjEquivDirectSecondPart->GetPDGMass();
  }
  G4double companionTotalEnergy =
    companionM0 + projectileKinEnergy - adjointPrimKinEnergy;
  G4double companionP2 =
    companionTotalEnergy * companionTotalEnergy - companionM0 * companionM0;

  // Projectile momentum
  G4double P_parallel =
    (adjointPrimP * adjointPrimP + projectileP2 - companionP2) /
    (2. * adjointPrimP);
  G4double P_perp = std::sqrt(projectileP2 - P_parallel * P_parallel);
  G4ThreeVector dir_parallel = theAdjointPrimary->GetMomentumDirection();
  G4double phi               = G4UniformRand() * twopi;
  G4ThreeVector projectileMomentum =
    G4ThreeVector(P_perp * std::cos(phi), P_perp * std::sin(phi), P_parallel);
  projectileMomentum.rotateUz(dir_parallel);

  if(!isScatProjToProj)
  {  // kill the primary and add a secondary
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    fParticleChange->AddSecondary(
      new G4DynamicParticle(fAdjEquivDirectPrimPart, projectileMomentum));
  }
  else
  {
    fParticleChange->ProposeEnergy(projectileKinEnergy);
    fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
  }
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointIonIonisationModel::DiffCrossSectionPerAtomPrimToSecond(
  G4double kinEnergyProj, G4double kinEnergyProd, G4double Z, G4double A)
{
  // Probably that here the Bragg Model should be also used for
  // kinEnergyProj/nuc<2MeV
  G4double dSigmadEprod = 0.;
  G4double Emax_proj    = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double Emin_proj    = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);

  G4double kinEnergyProjScaled = fMassRatio * kinEnergyProj;

  // the produced particle should have a kinetic energy smaller than the
  // projectile
  if(kinEnergyProj > Emin_proj && kinEnergyProj <= Emax_proj)
  {
    G4double Tmax = kinEnergyProj;

    G4double E1 = kinEnergyProd;
    G4double E2 = kinEnergyProd * 1.000001;
    G4double dE = (E2 - E1);
    G4double sigma1, sigma2;
    fDirectModel = fBraggIonDirectEMModel;
    if(kinEnergyProjScaled > 2. * MeV && !fUseOnlyBragg)
      fDirectModel = fBetheBlochDirectEMModel;
    sigma1 = fDirectModel->ComputeCrossSectionPerAtom(
      fDirectPrimaryPart, kinEnergyProj, Z, A, E1, 1.e20);
    sigma2 = fDirectModel->ComputeCrossSectionPerAtom(
      fDirectPrimaryPart, kinEnergyProj, Z, A, E2, 1.e20);

    dSigmadEprod = (sigma1 - sigma2) / dE;

    if(dSigmadEprod > 1.)
    {
      G4cout << "sigma1 " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << sigma1 << G4endl;
      G4cout << "sigma2 " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << sigma2 << G4endl;
      G4cout << "dsigma " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << dSigmadEprod << G4endl;
    }

    if(fDirectModel == fBetheBlochDirectEMModel)
    {
      // correction of differential cross section at high energy to correct for
      // the suppression of particle at secondary at high energy used in the
      // Bethe Bloch Model. This correction consist to multiply by g the
      // probability function used to test the rejection of a secondary Source
      // code taken from G4BetheBlochModel::SampleSecondaries

      G4double deltaKinEnergy = kinEnergyProd;

      G4double x = fFormFact * deltaKinEnergy;
      if(x > 1.e-6)
      {
        G4double totEnergy = kinEnergyProj + fMass;
        G4double etot2     = totEnergy * totEnergy;
        G4double beta2 = kinEnergyProj * (kinEnergyProj + 2.0 * fMass) / etot2;
        G4double f1    = 0.0;
        G4double f     = 1.0 - beta2 * deltaKinEnergy / Tmax;
        if(0.5 == fSpin)
        {
          f1 = 0.5 * deltaKinEnergy * deltaKinEnergy / etot2;
          f += f1;
        }
        G4double x1 = 1.0 + x;
        G4double gg = 1.0 / (x1 * x1);
        if(0.5 == fSpin)
        {
          G4double x2 =
            0.5 * electron_mass_c2 * deltaKinEnergy / (fMass * fMass);
          gg *= (1.0 + fMagMoment2 * (x2 - f1 / f) / (1.0 + x2));
        }
        if(gg > 1.0)
        {
          G4cout << "### G4BetheBlochModel in Adjoint Sim WARNING: gg= " << gg
                 << G4endl;
          gg = 1.;
        }
        dSigmadEprod *= gg;
      }
    }
  }
  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointIonIonisationModel::SetIon(G4ParticleDefinition* adj_ion,
                                         G4ParticleDefinition* fwd_ion)
{
  fDirectPrimaryPart      = fwd_ion;
  fAdjEquivDirectPrimPart = adj_ion;

  DefineProjectileProperty();
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointIonIonisationModel::CorrectPostStepWeight(
  G4ParticleChange* fParticleChange, G4double old_weight,
  G4double adjointPrimKinEnergy, G4double projectileKinEnergy, G4bool)
{
  // It is needed because the direct cross section used to compute the
  // differential cross section is not the one used in
  // the direct model where the GenericIon stuff is considered with correction
  // of effective charge.  In the direct model the samnepl of secondaries does
  // not reflect the integral cross section. The integral fwd cross section that
  // we used to compute the differential CS match the sample of secondaries in
  // the forward case despite the fact that its is not the same total CS than in
  // the FWD case. For this reason an extra weight correction is needed at the
  // end.

  G4double new_weight = old_weight;

  // the correction of CS due to the problem explained above
  G4double kinEnergyProjScaled = fMassRatio * projectileKinEnergy;
  fDirectModel                 = fBraggIonDirectEMModel;
  if(kinEnergyProjScaled > 2. * MeV && !fUseOnlyBragg)
    fDirectModel = fBetheBlochDirectEMModel;
  G4double UsedFwdCS = fDirectModel->ComputeCrossSectionPerAtom(
    fDirectPrimaryPart, projectileKinEnergy, 1, 1, fTcutSecond, 1.e20);
  G4double chargeSqRatio = 1.;
  if(fChargeSquare > 1.)
    chargeSqRatio = fDirectModel->GetChargeSquareRatio(
      fDirectPrimaryPart, fCurrentMaterial, projectileKinEnergy);
  G4double CorrectFwdCS =
    chargeSqRatio * fDirectModel->ComputeCrossSectionPerAtom(
                      G4GenericIon::GenericIon(), kinEnergyProjScaled, 1, 1,
                      fTcutSecond, 1.e20);
  // May be some check is needed if UsedFwdCS ==0 probably that then we should
  // avoid a secondary to be produced,
  if(UsedFwdCS > 0.)
    new_weight *= CorrectFwdCS / UsedFwdCS;

  // additional CS correction needed for cross section biasing in general.
  // May be wrong for ions. Most of the time not used.
  new_weight *=
    G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection() /
    fCsBiasingFactor;

  new_weight *= projectileKinEnergy / adjointPrimKinEnergy;

  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointIonIonisationModel::DefineProjectileProperty()
{
  // Slightly modified code taken from G4BetheBlochModel::SetParticle
  G4String pname = fDirectPrimaryPart->GetParticleName();

  fMass           = fDirectPrimaryPart->GetPDGMass();
  fMassRatio      = G4GenericIon::GenericIon()->GetPDGMass() / fMass;
  fSpin           = fDirectPrimaryPart->GetPDGSpin();
  G4double q      = fDirectPrimaryPart->GetPDGCharge() / eplus;
  fChargeSquare   = q * q;
  fRatio          = electron_mass_c2 / fMass;
  fOnePlusRatio2  = (1. + fRatio) * (1. + fRatio);
  fOneMinusRatio2 = (1. - fRatio) * (1. - fRatio);
  G4double magmom = fDirectPrimaryPart->GetPDGMagneticMoment() * fMass /
                    (0.5 * eplus * hbar_Planck * c_squared);
  fMagMoment2 = magmom * magmom - 1.0;
  if(fDirectPrimaryPart->GetLeptonNumber() == 0)
  {
    G4double x = 0.8426 * GeV;
    if(fSpin == 0.0 && fMass < GeV)
    {
      x = 0.736 * GeV;
    }
    else if(fMass > GeV)
    {
      x /= G4NistManager::Instance()->GetZ13(fMass / proton_mass_c2);
    }
    fFormFact = 2.0 * electron_mass_c2 / (x * x);
  }
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMaxForScatProjToProj(
  G4double primAdjEnergy)
{
  return primAdjEnergy * fOnePlusRatio2 /
         (fOneMinusRatio2 - 2. * fRatio * primAdjEnergy / fMass);
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMinForScatProjToProj(
  G4double primAdjEnergy, G4double tcut)
{
  return primAdjEnergy + tcut;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMaxForProdToProj(
  G4double)
{
  return GetHighEnergyLimit();
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMinForProdToProj(
  G4double primAdjEnergy)
{
  return (2. * primAdjEnergy - 4. * fMass +
          std::sqrt(4. * primAdjEnergy * primAdjEnergy + 16. * fMass * fMass +
                    8. * primAdjEnergy * fMass * (1. / fRatio + fRatio))) /
         4.;
}
