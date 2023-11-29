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

#include "G4AdjointBremsstrahlungModel.hh"

#include "G4AdjointCSManager.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Electron.hh"
#include "G4EmModelManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleChange.hh"
#include "G4PhysicalConstants.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackStatus.hh"

////////////////////////////////////////////////////////////////////////////////
G4AdjointBremsstrahlungModel::G4AdjointBremsstrahlungModel(G4VEmModel* aModel)
  : G4VEmAdjointModel("AdjointeBremModel")
{
  fDirectModel = aModel;
  Initialize();
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointBremsstrahlungModel::G4AdjointBremsstrahlungModel()
  : G4VEmAdjointModel("AdjointeBremModel")
{
  fDirectModel = new G4SeltzerBergerModel();
  Initialize();
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointBremsstrahlungModel::Initialize()
{
  SetUseMatrix(false);
  SetUseMatrixPerElement(false);

  fEmModelManagerForFwdModels = new G4EmModelManager();
  fEmModelManagerForFwdModels->AddEmModel(1, fDirectModel, nullptr, nullptr);
  SetApplyCutInRange(true);

  fElectron = G4Electron::Electron();
  fGamma    = G4Gamma::Gamma();

  fAdjEquivDirectPrimPart   = G4AdjointElectron::AdjointElectron();
  fAdjEquivDirectSecondPart = G4AdjointGamma::AdjointGamma();
  fDirectPrimaryPart        = fElectron;
  fSecondPartSameType       = false;

  fCSManager = G4AdjointCSManager::GetAdjointCSManager();
}

////////////////////////////////////////////////////////////////////////////////
G4AdjointBremsstrahlungModel::~G4AdjointBremsstrahlungModel()
{
  if(fEmModelManagerForFwdModels)
    delete fEmModelManagerForFwdModels;
}

////////////////////////////////////////////////////////////////////////////////
void G4AdjointBremsstrahlungModel::SampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  if(!fUseMatrix)
    return RapidSampleSecondaries(aTrack, isScatProjToProj, fParticleChange);

  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();
  DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());

  G4double adjointPrimKinEnergy   = theAdjointPrimary->GetKineticEnergy();
  G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();

  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  G4double projectileKinEnergy =
    SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy, isScatProjToProj);

  // Weight correction
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(),
                        adjointPrimKinEnergy, projectileKinEnergy,
                        isScatProjToProj);

  // Kinematic
  G4double projectileM0          = fAdjEquivDirectPrimPart->GetPDGMass();
  G4double projectileTotalEnergy = projectileM0 + projectileKinEnergy;
  G4double projectileP2 =
    projectileTotalEnergy * projectileTotalEnergy - projectileM0 * projectileM0;
  G4double projectileP = std::sqrt(projectileP2);

  // Angle of the gamma direction with the projectile taken from
  // G4eBremsstrahlungModel
  G4double u;
  if(0.25 > G4UniformRand())
    u = -std::log(G4UniformRand() * G4UniformRand()) / 0.625;
  else
    u = -std::log(G4UniformRand() * G4UniformRand()) / 1.875;

  G4double theta = u * electron_mass_c2 / projectileTotalEnergy;
  G4double sint  = std::sin(theta);
  G4double cost  = std::cos(theta);

  G4double phi = twopi * G4UniformRand();

  G4ThreeVector projectileMomentum =
    G4ThreeVector(std::cos(phi) * sint, std::sin(phi) * sint, cost) *
    projectileP;  // gamma frame
  if(isScatProjToProj)
  {  // the adjoint primary is the scattered e-
    G4ThreeVector gammaMomentum =
      (projectileTotalEnergy - adjointPrimTotalEnergy) *
      G4ThreeVector(0., 0., 1.);
    G4ThreeVector dirProd = projectileMomentum - gammaMomentum;
    G4double cost1        = std::cos(dirProd.angle(projectileMomentum));
    G4double sint1        = std::sqrt(1. - cost1 * cost1);
    projectileMomentum =
      G4ThreeVector(std::cos(phi) * sint1, std::sin(phi) * sint1, cost1) *
      projectileP;
  }

  projectileMomentum.rotateUz(theAdjointPrimary->GetMomentumDirection());

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
void G4AdjointBremsstrahlungModel::RapidSampleSecondaries(
  const G4Track& aTrack, G4bool isScatProjToProj,
  G4ParticleChange* fParticleChange)
{
  const G4DynamicParticle* theAdjointPrimary = aTrack.GetDynamicParticle();
  DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());

  G4double adjointPrimKinEnergy   = theAdjointPrimary->GetKineticEnergy();
  G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();

  if(adjointPrimKinEnergy > GetHighEnergyLimit() * 0.999)
  {
    return;
  }

  G4double projectileKinEnergy = 0.;
  G4double gammaEnergy         = 0.;
  G4double diffCSUsed          = 0.;
  if(!isScatProjToProj)
  {
    gammaEnergy   = adjointPrimKinEnergy;
    G4double Emax = GetSecondAdjEnergyMaxForProdToProj(adjointPrimKinEnergy);
    G4double Emin = GetSecondAdjEnergyMinForProdToProj(adjointPrimKinEnergy);
    if(Emin >= Emax)
      return;
    projectileKinEnergy = Emin * std::pow(Emax / Emin, G4UniformRand());
    diffCSUsed          = fCsBiasingFactor * fLastCZ / projectileKinEnergy;
  }
  else
  {
    G4double Emax =
      GetSecondAdjEnergyMaxForScatProjToProj(adjointPrimKinEnergy);
    G4double Emin =
      GetSecondAdjEnergyMinForScatProjToProj(adjointPrimKinEnergy, fTcutSecond);
    if(Emin >= Emax)
      return;
    G4double f1 = (Emin - adjointPrimKinEnergy) / Emin;
    G4double f2 = (Emax - adjointPrimKinEnergy) / Emax / f1;
    projectileKinEnergy =
      adjointPrimKinEnergy / (1. - f1 * std::pow(f2, G4UniformRand()));
    gammaEnergy = projectileKinEnergy - adjointPrimKinEnergy;
    diffCSUsed =
      fLastCZ * adjointPrimKinEnergy / projectileKinEnergy / gammaEnergy;
  }

  // Weight correction:
  // First w_corr is set to the ratio between adjoint total CS and fwd total CS
  // if this has to be done in the model.
  // For the case of forced interaction this will be done in the PostStepDoIt of
  // the forced interaction.  It is important to set the weight before the
  // creation of the secondary
  G4double w_corr = fOutsideWeightFactor;
  if(fInModelWeightCorr)
  {
    w_corr = fCSManager->GetPostStepWeightCorrection();
  }

  // Then another correction is needed due to the fact that a biaised
  // differential CS has been used rather than the one consistent with the
  // direct model Here we consider the true diffCS as the one obtained by the
  // numerical differentiation over Tcut of the direct CS, corrected by the
  // Migdal term. Basically any other differential CS could be used here
  // (example Penelope).
  G4double diffCS = DiffCrossSectionPerVolumePrimToSecond(
    fCurrentMaterial, projectileKinEnergy, gammaEnergy);
  w_corr *= diffCS / diffCSUsed;

  G4double new_weight = aTrack.GetWeight() * w_corr;
  fParticleChange->SetParentWeightByProcess(false);
  fParticleChange->SetSecondaryWeightByProcess(false);
  fParticleChange->ProposeParentWeight(new_weight);

  // Kinematic
  G4double projectileM0          = fAdjEquivDirectPrimPart->GetPDGMass();
  G4double projectileTotalEnergy = projectileM0 + projectileKinEnergy;
  G4double projectileP2 =
    projectileTotalEnergy * projectileTotalEnergy - projectileM0 * projectileM0;
  G4double projectileP = std::sqrt(projectileP2);

  // Use the angular model of the forward model to generate the gamma direction
  // Dummy dynamic particle to use the model
  G4DynamicParticle* aDynPart =
    new G4DynamicParticle(fElectron, G4ThreeVector(0., 0., 1.) * projectileP);

  // Get the element from the direct model
  const G4Element* elm = fDirectModel->SelectRandomAtom(
    fCurrentCouple, fElectron, projectileKinEnergy, fTcutSecond);
  G4int Z         = elm->GetZasInt();
  G4double energy = aDynPart->GetTotalEnergy() - gammaEnergy;
  G4ThreeVector projectileMomentum =
    fDirectModel->GetAngularDistribution()->SampleDirection(aDynPart, energy, Z,
                                            fCurrentMaterial) * projectileP;
  G4double phi = projectileMomentum.getPhi();

  if(isScatProjToProj)
  {  // the adjoint primary is the scattered e-
    G4ThreeVector gammaMomentum =
      (projectileTotalEnergy - adjointPrimTotalEnergy) *
      G4ThreeVector(0., 0., 1.);
    G4ThreeVector dirProd = projectileMomentum - gammaMomentum;
    G4double cost1        = std::cos(dirProd.angle(projectileMomentum));
    G4double sint1        = std::sqrt(1. - cost1 * cost1);
    projectileMomentum =
      G4ThreeVector(std::cos(phi) * sint1, std::sin(phi) * sint1, cost1) *
      projectileP;
  }

  projectileMomentum.rotateUz(theAdjointPrimary->GetMomentumDirection());

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
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond(
  const G4Material* aMaterial,
  G4double kinEnergyProj,  // kin energy of primary before interaction
  G4double kinEnergyProd   // kinetic energy of the secondary particle
)
{
  if(!fIsDirectModelInitialised)
  {
    fEmModelManagerForFwdModels->Initialise(fElectron, fGamma, 0);
    fIsDirectModelInitialised = true;
  }
  return G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToSecond(
    aMaterial, kinEnergyProj, kinEnergyProd);
}

////////////////////////////////////////////////////////////////////////////////
G4double G4AdjointBremsstrahlungModel::AdjointCrossSection(
  const G4MaterialCutsCouple* aCouple, G4double primEnergy,
  G4bool isScatProjToProj)
{
  static constexpr G4double maxEnergy = 100. * MeV / 2.718281828459045;
  // 2.78.. == std::exp(1.)
  if(!fIsDirectModelInitialised)
  {
    fEmModelManagerForFwdModels->Initialise(fElectron, fGamma, 0);
    fIsDirectModelInitialised = true;
  }
  if(fUseMatrix)
    return G4VEmAdjointModel::AdjointCrossSection(aCouple, primEnergy,
                                                  isScatProjToProj);
  DefineCurrentMaterial(aCouple);
  G4double Cross = 0.;
  // this gives the constant above
  fLastCZ = fDirectModel->CrossSectionPerVolume(
    aCouple->GetMaterial(), fDirectPrimaryPart, 100. * MeV, maxEnergy);

  if(!isScatProjToProj)
  {
    G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProj(primEnergy);
    G4double Emin_proj = GetSecondAdjEnergyMinForProdToProj(primEnergy);
    if(Emax_proj > Emin_proj && primEnergy > fTcutSecond)
      Cross = fCsBiasingFactor * fLastCZ * std::log(Emax_proj / Emin_proj);
  }
  else
  {
    G4double Emax_proj = GetSecondAdjEnergyMaxForScatProjToProj(primEnergy);
    G4double Emin_proj =
      GetSecondAdjEnergyMinForScatProjToProj(primEnergy, fTcutSecond);
    if(Emax_proj > Emin_proj)
      Cross = fLastCZ * std::log((Emax_proj - primEnergy) * Emin_proj /
                                 Emax_proj / (Emin_proj - primEnergy));
  }
  return Cross;
}
