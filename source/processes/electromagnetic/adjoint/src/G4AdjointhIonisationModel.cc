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

#include "G4AdjointhIonisationModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointProton.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4NistManager.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4Proton.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"

////////////////////////////////////////////////////////////////////////////////
G4AdjointhIonisationModel::G4AdjointhIonisationModel(G4ParticleDefinition* pDef)
  : G4VEmAdjointModel("Adjoint_hIonisation")
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

  fDirectModel              = new G4BetheBlochModel(pDef);
  fBraggDirectEMModel       = new G4BraggModel(pDef);
  fAdjEquivDirectSecondPart = G4AdjointElectron::AdjointElectron();
  fDirectPrimaryPart        = pDef;

  if(pDef == G4Proton::Proton())
  {
    fAdjEquivDirectPrimPart = G4AdjointProton::AdjointProton();
  }

  DefineProjectileProperty();
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointhIonisationModel::~G4AdjointhIonisationModel() {}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointhIonisationModel::SampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  if(!fUseMatrix)
    return RapidSampleSecondaries(aTrack, isScatProjToProj, fParticleChange);

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
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(),
                        adjointPrimKinEnergy, projectileKinEnergy,
                        isScatProjToProj);
  // Caution!!! this weight correction should be always applied

  // Kinematic:
  // we consider a two body elastic scattering for the forward processes where
  // the projectile knock on an e- at rest and gives it part of its energy
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
void G4AdjointhIonisationModel::RapidSampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();
  DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());

  G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
  G4double adjointPrimP         = theAdjointPrimary->GetTotalMomentum();

  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  G4double projectileKinEnergy = 0.;
  G4double eEnergy             = 0.;
  G4double newCS =
    fCurrentMaterial->GetElectronDensity() * twopi_mc2_rcl2 * fMass;
  if(!isScatProjToProj)
  {  // 1/E^2 distribution

    eEnergy       = adjointPrimKinEnergy;
    G4double Emax = GetSecondAdjEnergyMaxForProdToProj(adjointPrimKinEnergy);
    G4double Emin = GetSecondAdjEnergyMinForProdToProj(adjointPrimKinEnergy);
    if(Emin >= Emax)
      return;
    G4double a = 1. / Emax;
    G4double b = 1. / Emin;
    newCS      = newCS * (b - a) / eEnergy;

    projectileKinEnergy = 1. / (b - (b - a) * G4UniformRand());
  }
  else
  {
    G4double Emax =
      GetSecondAdjEnergyMaxForScatProjToProj(adjointPrimKinEnergy);
    G4double Emin =
      GetSecondAdjEnergyMinForScatProjToProj(adjointPrimKinEnergy, fTcutSecond);
    if(Emin >= Emax)
      return;
    G4double diff1 = Emin - adjointPrimKinEnergy;
    G4double diff2 = Emax - adjointPrimKinEnergy;

    G4double t1    = adjointPrimKinEnergy * (1. / diff1 - 1. / diff2);
    G4double t2    = adjointPrimKinEnergy * (1. / Emin - 1. / Emax);
    G4double t3    = 2. * std::log(Emax / Emin);
    G4double sum_t = t1 + t2 + t3;
    newCS      = newCS * sum_t / adjointPrimKinEnergy / adjointPrimKinEnergy;
    G4double t = G4UniformRand() * sum_t;
    if(t <= t1)
    {
      G4double q          = G4UniformRand() * t1 / adjointPrimKinEnergy;
      projectileKinEnergy = adjointPrimKinEnergy + 1. / (1. / diff1 - q);
    }
    else if(t <= t2)
    {
      G4double q          = G4UniformRand() * t2 / adjointPrimKinEnergy;
      projectileKinEnergy = 1. / (1. / Emin - q);
    }
    else
    {
      projectileKinEnergy = Emin * std::pow(Emax / Emin, G4UniformRand());
    }
    eEnergy = projectileKinEnergy - adjointPrimKinEnergy;
  }

  G4double diffCS_perAtom_Used = twopi_mc2_rcl2 * fMass * adjointPrimKinEnergy /
                                 projectileKinEnergy / projectileKinEnergy /
                                 eEnergy / eEnergy;

  // Weight correction
  // First w_corr is set to the ratio between adjoint total CS and fwd total CS
  G4double w_corr =
    G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();

  w_corr *= newCS / fLastCS;
  // Then another correction is needed due to the fact that a biaised
  // differential CS has been used rather than the one consistent with the
  // direct model. Here we consider the true diffCS as the one obtained by the
  // numerical differentiation over Tcut of the direct CS

  G4double diffCS =
    DiffCrossSectionPerAtomPrimToSecond(projectileKinEnergy, eEnergy, 1, 1);
  w_corr *= diffCS / diffCS_perAtom_Used;

  if (isScatProjToProj && fTcutSecond>0.005) w_corr=1.;

  G4double new_weight = aTrack.GetWeight() * w_corr;
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);

  // Kinematic:
  // we consider a two body elastic scattering for the forward processes where
  // the projectile knocks on an e- at rest and gives it part of its energy
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
G4double G4AdjointhIonisationModel::DiffCrossSectionPerAtomPrimToSecond(
  G4double kinEnergyProj, G4double kinEnergyProd, G4double Z, G4double A)
{  // Probably here the Bragg Model should be also used for
   // kinEnergyProj/nuc < 2 MeV

  G4double dSigmadEprod = 0.;
  G4double Emax_proj    = GetSecondAdjEnergyMaxForProdToProj(kinEnergyProd);
  G4double Emin_proj    = GetSecondAdjEnergyMinForProdToProj(kinEnergyProd);

  // the produced particle should have a kinetic energy smaller than the
  // projectile
  if(kinEnergyProj > Emin_proj && kinEnergyProj <= Emax_proj)
  {
    G4double Tmax = kinEnergyProj;
    G4double E1   = kinEnergyProd;
    //1.0006 factor seems to give the best diff CS, important impact on proton correction factor
    G4double E2   = kinEnergyProd *1.0006;
    G4double sigma1, sigma2;
    if(kinEnergyProj > 2. * MeV)
    {
      sigma1 = fDirectModel->ComputeCrossSectionPerAtom(
        fDirectPrimaryPart, kinEnergyProj, Z, A, E1, 1.e20);
      sigma2 = fDirectModel->ComputeCrossSectionPerAtom(
        fDirectPrimaryPart, kinEnergyProj, Z, A, E2, 1.e20);
    }
    else
    {
      sigma1 = fBraggDirectEMModel->ComputeCrossSectionPerAtom(
        fDirectPrimaryPart, kinEnergyProj, Z, A, E1, 1.e20);
      sigma2 = fBraggDirectEMModel->ComputeCrossSectionPerAtom(
        fDirectPrimaryPart, kinEnergyProj, Z, A, E2, 1.e20);
    }

    dSigmadEprod = (sigma1 - sigma2) / (E2 - E1);
    if(dSigmadEprod > 1.)
    {
      G4cout << "sigma1 " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << sigma1 << G4endl;
      G4cout << "sigma2 " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << sigma2 << G4endl;
      G4cout << "dsigma " << kinEnergyProj / MeV << '\t' << kinEnergyProd / MeV
             << '\t' << dSigmadEprod << G4endl;
    }

    // correction of differential cross section at high energy to correct for
    // the suppression of particle at secondary at high energy used in the Bethe
    // Bloch Model. This correction consists of multiplying by g the probability
    // function used to test the rejection of a secondary. Source code taken
    // from G4BetheBlochModel::SampleSecondaries
    G4double deltaKinEnergy = kinEnergyProd;

    // projectile formfactor - suppression of high energy
    // delta-electron production at high energy
    G4double x = fFormFact * deltaKinEnergy;
    if(x > 1.e-6)
    {
      G4double totEnergy = kinEnergyProj + fMass;
      G4double etot2     = totEnergy * totEnergy;
      G4double beta2 = kinEnergyProj * (kinEnergyProj + 2.0 * fMass) / etot2;
      G4double f     = 1.0 - beta2 * deltaKinEnergy / Tmax;
      G4double f1    = 0.0;
      if(0.5 == fSpin)
      {
        f1 = 0.5 * deltaKinEnergy * deltaKinEnergy / etot2;
        f += f1;
      }
      G4double x1 = 1.0 + x;
      G4double gg = 1.0 / (x1 * x1);
      if(0.5 == fSpin)
      {
        G4double x2 = 0.5 * electron_mass_c2 * deltaKinEnergy / (fMass * fMass);
        gg *= (1.0 + fMagMoment2 * (x2 - f1 / f) / (1.0 + x2));
      }
      if(gg > 1.0)
      {
        G4cout << "### G4BetheBlochModel in Adjoint Sim WARNING: g= " << g
               << G4endl;
        gg = 1.;
      }
      dSigmadEprod *= gg;
    }
  }

  return dSigmadEprod;
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointhIonisationModel::DefineProjectileProperty()
{
  // Slightly modified code taken from G4BetheBlochModel::SetParticle
  G4String pname = fDirectPrimaryPart->GetParticleName();

  fMass           = fDirectPrimaryPart->GetPDGMass();
  fSpin           = fDirectPrimaryPart->GetPDGSpin();
  fMassRatio      = electron_mass_c2 / fMass;
  fOnePlusRatio2  = (1. + fMassRatio) * (1. + fMassRatio);
  fOneMinusRatio2 = (1. - fMassRatio) * (1. - fMassRatio);
  G4double magmom = fDirectPrimaryPart->GetPDGMagneticMoment() * fMass /
                    (0.5 * eplus * hbar_Planck * c_squared);
  fMagMoment2 = magmom * magmom - 1.0;
  fFormFact   = 0.0;
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

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointhIonisationModel::AdjointCrossSection(
  const G4MaterialCutsCouple* aCouple, G4double primEnergy,
  G4bool isScatProjToProj)
{
  if(fUseMatrix)
    return G4VEmAdjointModel::AdjointCrossSection(aCouple, primEnergy,
                                                  isScatProjToProj);
  DefineCurrentMaterial(aCouple);

  G4double Cross =
    fCurrentMaterial->GetElectronDensity() * twopi_mc2_rcl2 * fMass;

  if(!isScatProjToProj)
  {
    G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProj(primEnergy);
    G4double Emin_proj = GetSecondAdjEnergyMinForProdToProj(primEnergy);
    if(Emax_proj > Emin_proj && primEnergy > fTcutSecond)
    {
      Cross *= (1. / Emin_proj - 1. / Emax_proj) / primEnergy;
    }
    else
      Cross = 0.;
  }
  else
  {
    G4double Emax_proj = GetSecondAdjEnergyMaxForScatProjToProj(primEnergy);
    G4double Emin_proj =
      GetSecondAdjEnergyMinForScatProjToProj(primEnergy, fTcutSecond);
    G4double diff1 = Emin_proj - primEnergy;
    G4double diff2 = Emax_proj - primEnergy;
    G4double t1 =
      (1. / diff1 + 1. / Emin_proj - 1. / diff2 - 1. / Emax_proj) / primEnergy;
    G4double t2 =
      2. * std::log(Emax_proj / Emin_proj) / primEnergy / primEnergy;
    Cross *= (t1 + t2);
  }
  fLastCS = Cross;
  return Cross;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMaxForScatProjToProj(
  G4double primAdjEnergy)
{
  G4double Tmax = primAdjEnergy * fOnePlusRatio2 /
                  (fOneMinusRatio2 - 2. * fMassRatio * primAdjEnergy / fMass);
  return Tmax;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMinForScatProjToProj(
  G4double primAdjEnergy, G4double tcut)
{
  return primAdjEnergy + tcut;
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMaxForProdToProj(G4double)
{
  return GetHighEnergyLimit();
}

//////////////////////////////////////////////////////////////////////////////
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMinForProdToProj(
  G4double primAdjEnergy)
{
  G4double Tmin =
    (2. * primAdjEnergy - 4. * fMass +
     std::sqrt(4. * primAdjEnergy * primAdjEnergy + 16. * fMass * fMass +
               8. * primAdjEnergy * fMass * (1. / fMassRatio + fMassRatio))) /
    4.;
  return Tmin;
}
