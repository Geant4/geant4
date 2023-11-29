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
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation,
//      21-5-98 V.Grichine
//      28-05-01, V.Ivanchenko minor changes to provide ANSI -wall compilation
//      04.03.05, V.Grichine: get local field interface
//      19-05-06, V.Ivanchenko rename from G4SynchrotronRadiation
//
///////////////////////////////////////////////////////////////////////////

#include "G4SynchrotronRadiationInMat.hh"

#include "G4EmProcessSubType.hh"
#include "G4Field.hh"
#include "G4FieldManager.hh"
#include "G4Integrator.hh"
#include "G4PhysicalConstants.hh"
#include "G4PropagatorInField.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsModelCatalog.hh"

const G4double G4SynchrotronRadiationInMat::fIntegralProbabilityOfSR[200] = {
  1.000000e+00, 9.428859e-01, 9.094095e-01, 8.813971e-01, 8.565154e-01,
  8.337008e-01, 8.124961e-01, 7.925217e-01, 7.735517e-01, 7.554561e-01,
  7.381233e-01, 7.214521e-01, 7.053634e-01, 6.898006e-01, 6.747219e-01,
  6.600922e-01, 6.458793e-01, 6.320533e-01, 6.185872e-01, 6.054579e-01,
  5.926459e-01, 5.801347e-01, 5.679103e-01, 5.559604e-01, 5.442736e-01,
  5.328395e-01, 5.216482e-01, 5.106904e-01, 4.999575e-01, 4.894415e-01,
  4.791351e-01, 4.690316e-01, 4.591249e-01, 4.494094e-01, 4.398800e-01,
  4.305320e-01, 4.213608e-01, 4.123623e-01, 4.035325e-01, 3.948676e-01,
  3.863639e-01, 3.780179e-01, 3.698262e-01, 3.617858e-01, 3.538933e-01,
  3.461460e-01, 3.385411e-01, 3.310757e-01, 3.237474e-01, 3.165536e-01,
  3.094921e-01, 3.025605e-01, 2.957566e-01, 2.890784e-01, 2.825237e-01,
  2.760907e-01, 2.697773e-01, 2.635817e-01, 2.575020e-01, 2.515365e-01,
  2.456834e-01, 2.399409e-01, 2.343074e-01, 2.287812e-01, 2.233607e-01,
  2.180442e-01, 2.128303e-01, 2.077174e-01, 2.027040e-01, 1.977885e-01,
  1.929696e-01, 1.882457e-01, 1.836155e-01, 1.790775e-01, 1.746305e-01,
  1.702730e-01, 1.660036e-01, 1.618212e-01, 1.577243e-01, 1.537117e-01,
  1.497822e-01, 1.459344e-01, 1.421671e-01, 1.384791e-01, 1.348691e-01,
  1.313360e-01, 1.278785e-01, 1.244956e-01, 1.211859e-01, 1.179483e-01,
  1.147818e-01, 1.116850e-01, 1.086570e-01, 1.056966e-01, 1.028026e-01,
  9.997405e-02, 9.720975e-02, 9.450865e-02, 9.186969e-02, 8.929179e-02,
  8.677391e-02, 8.431501e-02, 8.191406e-02, 7.957003e-02, 7.728192e-02,
  7.504872e-02, 7.286944e-02, 7.074311e-02, 6.866874e-02, 6.664538e-02,
  6.467208e-02, 6.274790e-02, 6.087191e-02, 5.904317e-02, 5.726079e-02,
  5.552387e-02, 5.383150e-02, 5.218282e-02, 5.057695e-02, 4.901302e-02,
  4.749020e-02, 4.600763e-02, 4.456450e-02, 4.315997e-02, 4.179325e-02,
  4.046353e-02, 3.917002e-02, 3.791195e-02, 3.668855e-02, 3.549906e-02,
  3.434274e-02, 3.321884e-02, 3.212665e-02, 3.106544e-02, 3.003452e-02,
  2.903319e-02, 2.806076e-02, 2.711656e-02, 2.619993e-02, 2.531021e-02,
  2.444677e-02, 2.360897e-02, 2.279620e-02, 2.200783e-02, 2.124327e-02,
  2.050194e-02, 1.978324e-02, 1.908662e-02, 1.841151e-02, 1.775735e-02,
  1.712363e-02, 1.650979e-02, 1.591533e-02, 1.533973e-02, 1.478250e-02,
  1.424314e-02, 1.372117e-02, 1.321613e-02, 1.272755e-02, 1.225498e-02,
  1.179798e-02, 1.135611e-02, 1.092896e-02, 1.051609e-02, 1.011712e-02,
  9.731635e-03, 9.359254e-03, 8.999595e-03, 8.652287e-03, 8.316967e-03,
  7.993280e-03, 7.680879e-03, 7.379426e-03, 7.088591e-03, 6.808051e-03,
  6.537491e-03, 6.276605e-03, 6.025092e-03, 5.782661e-03, 5.549027e-03,
  5.323912e-03, 5.107045e-03, 4.898164e-03, 4.697011e-03, 4.503336e-03,
  4.316896e-03, 4.137454e-03, 3.964780e-03, 3.798649e-03, 3.638843e-03,
  3.485150e-03, 3.337364e-03, 3.195284e-03, 3.058715e-03, 2.927469e-03,
  2.801361e-03, 2.680213e-03, 2.563852e-03, 2.452110e-03, 2.344824e-03
};

///////////////////////////////////////////////////////////////////////
//  Constructor
G4SynchrotronRadiationInMat::G4SynchrotronRadiationInMat(
  const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
  , theGamma(G4Gamma::Gamma())
  , theElectron(G4Electron::Electron())
  , thePositron(G4Positron::Positron())
  , LowestKineticEnergy(10. * keV)
  , fAlpha(0.0)
  , fRootNumber(80)
  , fVerboseLevel(verboseLevel)
{
  G4TransportationManager* transportMgr =
    G4TransportationManager::GetTransportationManager();

  fFieldPropagator = transportMgr->GetPropagatorInField();
  secID = G4PhysicsModelCatalog::GetModelID("model_SynchrotronRadiation");
  SetProcessSubType(fSynchrotronRadiation);
  CutInRange = GammaCutInKineticEnergyNow = ElectronCutInKineticEnergyNow =
    PositronCutInKineticEnergyNow = ParticleCutInKineticEnergyNow = fKsi =
      fPsiGamma = fEta = fOrderAngleK = 0.0;
}

/////////////////////////////////////////////////////////////////////////
// Destructor
G4SynchrotronRadiationInMat::~G4SynchrotronRadiationInMat() = default;

G4bool G4SynchrotronRadiationInMat::IsApplicable(
  const G4ParticleDefinition& particle)
{
  return ((&particle == (const G4ParticleDefinition*) theElectron) ||
          (&particle == (const G4ParticleDefinition*) thePositron));
}

G4double G4SynchrotronRadiationInMat::GetLambdaConst() { return fLambdaConst; }

G4double G4SynchrotronRadiationInMat::GetEnergyConst() { return fEnergyConst; }

// Production of synchrotron X-ray photon
// Geant4 internal units.
G4double G4SynchrotronRadiationInMat::GetMeanFreePath(
  const G4Track& trackData, G4double, G4ForceCondition* condition)
{
  // gives the MeanFreePath in GEANT4 internal units
  G4double MeanFreePath;

  const G4DynamicParticle* aDynamicParticle = trackData.GetDynamicParticle();

  *condition = NotForced;

  G4double gamma =
    aDynamicParticle->GetTotalEnergy() / aDynamicParticle->GetMass();

  G4double particleCharge = aDynamicParticle->GetDefinition()->GetPDGCharge();

  G4double KineticEnergy = aDynamicParticle->GetKineticEnergy();

  if(KineticEnergy < LowestKineticEnergy || gamma < 1.0e3)
    MeanFreePath = DBL_MAX;
  else
  {
    G4ThreeVector FieldValue;
    const G4Field* pField = nullptr;

    G4FieldManager* fieldMgr = nullptr;
    G4bool fieldExertsForce  = false;

    if((particleCharge != 0.0))
    {
      fieldMgr =
        fFieldPropagator->FindAndSetFieldManager(trackData.GetVolume());

      if(fieldMgr != nullptr)
      {
        // If the field manager has no field, there is no field !
        fieldExertsForce = (fieldMgr->GetDetectorField() != nullptr);
      }
    }
    if(fieldExertsForce)
    {
      pField                     = fieldMgr->GetDetectorField();
      G4ThreeVector globPosition = trackData.GetPosition();

      G4double globPosVec[4], FieldValueVec[6];

      globPosVec[0] = globPosition.x();
      globPosVec[1] = globPosition.y();
      globPosVec[2] = globPosition.z();
      globPosVec[3] = trackData.GetGlobalTime();

      pField->GetFieldValue(globPosVec, FieldValueVec);

      FieldValue =
        G4ThreeVector(FieldValueVec[0], FieldValueVec[1], FieldValueVec[2]);

      G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection();
      G4ThreeVector unitMcrossB  = FieldValue.cross(unitMomentum);
      G4double perpB             = unitMcrossB.mag();
      G4double beta              = aDynamicParticle->GetTotalMomentum() /
                      (aDynamicParticle->GetTotalEnergy());

      if(perpB > 0.0)
        MeanFreePath = fLambdaConst * beta / perpB;
      else
        MeanFreePath = DBL_MAX;
    }
    else
      MeanFreePath = DBL_MAX;
  }
  if(fVerboseLevel > 0)
  {
    G4cout << "G4SynchrotronRadiationInMat::MeanFreePath = " << MeanFreePath / m
           << " m" << G4endl;
  }
  return MeanFreePath;
}

////////////////////////////////////////////////////////////////////////////////
G4VParticleChange* G4SynchrotronRadiationInMat::PostStepDoIt(
  const G4Track& trackData, const G4Step& stepData)

{
  aParticleChange.Initialize(trackData);

  const G4DynamicParticle* aDynamicParticle = trackData.GetDynamicParticle();

  G4double gamma =
    aDynamicParticle->GetTotalEnergy() / (aDynamicParticle->GetMass());

  if(gamma <= 1.0e3)
  {
    return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
  }
  G4double particleCharge = aDynamicParticle->GetDefinition()->GetPDGCharge();

  G4ThreeVector FieldValue;
  const G4Field* pField = nullptr;

  G4FieldManager* fieldMgr = nullptr;
  G4bool fieldExertsForce  = false;

  if((particleCharge != 0.0))
  {
    fieldMgr = fFieldPropagator->FindAndSetFieldManager(trackData.GetVolume());
    if(fieldMgr != nullptr)
    {
      // If the field manager has no field, there is no field !
      fieldExertsForce = (fieldMgr->GetDetectorField() != nullptr);
    }
  }
  if(fieldExertsForce)
  {
    pField                     = fieldMgr->GetDetectorField();
    G4ThreeVector globPosition = trackData.GetPosition();
    G4double globPosVec[4], FieldValueVec[6];
    globPosVec[0] = globPosition.x();
    globPosVec[1] = globPosition.y();
    globPosVec[2] = globPosition.z();
    globPosVec[3] = trackData.GetGlobalTime();

    pField->GetFieldValue(globPosVec, FieldValueVec);
    FieldValue =
      G4ThreeVector(FieldValueVec[0], FieldValueVec[1], FieldValueVec[2]);

    G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection();
    G4ThreeVector unitMcrossB  = FieldValue.cross(unitMomentum);
    G4double perpB             = unitMcrossB.mag();
    if(perpB > 0.0)
    {
      // M-C of synchrotron photon energy
      G4double energyOfSR = GetRandomEnergySR(gamma, perpB);

      if(fVerboseLevel > 0)
      {
        G4cout << "SR photon energy = " << energyOfSR / keV << " keV" << G4endl;
      }

      // check against insufficient energy
      if(energyOfSR <= 0.0)
      {
        return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
      }
      G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();
      G4ParticleMomentum particleDirection =
        aDynamicParticle->GetMomentumDirection();

      // M-C of its direction, simplified dipole busted approach
      G4double cosTheta, sinTheta, fcos, beta;

      do
      {
        cosTheta = 1. - 2. * G4UniformRand();
        fcos     = (1 + cosTheta * cosTheta) * 0.5;
      }
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
      while(fcos < G4UniformRand());

      beta = std::sqrt(1. - 1. / (gamma * gamma));

      cosTheta = (cosTheta + beta) / (1. + beta * cosTheta);

      if(cosTheta > 1.)
        cosTheta = 1.;
      if(cosTheta < -1.)
        cosTheta = -1.;

      sinTheta = std::sqrt(1. - cosTheta * cosTheta);

      G4double Phi = twopi * G4UniformRand();

      G4double dirx = sinTheta * std::cos(Phi);
      G4double diry = sinTheta * std::sin(Phi);
      G4double dirz = cosTheta;

      G4ThreeVector gammaDirection(dirx, diry, dirz);
      gammaDirection.rotateUz(particleDirection);

      // polarization of new gamma
      G4ThreeVector gammaPolarization = FieldValue.cross(gammaDirection);
      gammaPolarization               = gammaPolarization.unit();

      // create G4DynamicParticle object for the SR photon
      auto aGamma =
        new G4DynamicParticle(G4Gamma::Gamma(), gammaDirection, energyOfSR);
      aGamma->SetPolarization(gammaPolarization.x(), gammaPolarization.y(),
                              gammaPolarization.z());

      aParticleChange.SetNumberOfSecondaries(1);

      // Update the incident particle
      G4double newKinEnergy = kineticEnergy - energyOfSR;

      if(newKinEnergy > 0.)
      {
        aParticleChange.ProposeMomentumDirection(particleDirection);
        aParticleChange.ProposeEnergy(newKinEnergy);
        aParticleChange.ProposeLocalEnergyDeposit(0.);
      }
      else
      {
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeLocalEnergyDeposit(0.);
        G4double charge = aDynamicParticle->GetDefinition()->GetPDGCharge();
        if(charge < 0.)
        {
          aParticleChange.ProposeTrackStatus(fStopAndKill);
        }
        else
        {
          aParticleChange.ProposeTrackStatus(fStopButAlive);
        }
      }
      
      // Create the G4Track
      G4Track* aSecondaryTrack = new G4Track(aGamma, trackData.GetGlobalTime(), trackData.GetPosition());
      aSecondaryTrack->SetTouchableHandle(stepData.GetPostStepPoint()->GetTouchableHandle());
      aSecondaryTrack->SetParentID(trackData.GetTrackID());
      aSecondaryTrack->SetCreatorModelID(secID);
      aParticleChange.AddSecondary(aSecondaryTrack);
    }
    else
    {
      return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(trackData, stepData);
}

G4double G4SynchrotronRadiationInMat::GetPhotonEnergy(const G4Track& trackData,
                                                      const G4Step&)
{
  G4int i;
  G4double energyOfSR = -1.0;

  const G4DynamicParticle* aDynamicParticle = trackData.GetDynamicParticle();

  G4double gamma =
    aDynamicParticle->GetTotalEnergy() / (aDynamicParticle->GetMass());

  G4double particleCharge = aDynamicParticle->GetDefinition()->GetPDGCharge();

  G4ThreeVector FieldValue;
  const G4Field* pField = nullptr;

  G4FieldManager* fieldMgr = nullptr;
  G4bool fieldExertsForce  = false;

  if((particleCharge != 0.0))
  {
    fieldMgr = fFieldPropagator->FindAndSetFieldManager(trackData.GetVolume());
    if(fieldMgr != nullptr)
    {
      // If the field manager has no field, there is no field !
      fieldExertsForce = (fieldMgr->GetDetectorField() != nullptr);
    }
  }
  if(fieldExertsForce)
  {
    pField                     = fieldMgr->GetDetectorField();
    G4ThreeVector globPosition = trackData.GetPosition();
    G4double globPosVec[3], FieldValueVec[3];

    globPosVec[0] = globPosition.x();
    globPosVec[1] = globPosition.y();
    globPosVec[2] = globPosition.z();

    pField->GetFieldValue(globPosVec, FieldValueVec);
    FieldValue =
      G4ThreeVector(FieldValueVec[0], FieldValueVec[1], FieldValueVec[2]);

    G4ThreeVector unitMomentum = aDynamicParticle->GetMomentumDirection();
    G4ThreeVector unitMcrossB  = FieldValue.cross(unitMomentum);
    G4double perpB             = unitMcrossB.mag();
    if(perpB > 0.0)
    {
      // M-C of synchrotron photon energy
      G4double random = G4UniformRand();
      for(i = 0; i < 200; ++i)
      {
        if(random >= fIntegralProbabilityOfSR[i])
          break;
      }
      energyOfSR = 0.0001 * i * i * fEnergyConst * gamma * gamma * perpB;

      // check against insufficient energy
      if(energyOfSR <= 0.0)
      {
        return -1.0;
      }
    }
    else
    {
      return -1.0;
    }
  }
  return energyOfSR;
}

/////////////////////////////////////////////////////////////////////////////////
G4double G4SynchrotronRadiationInMat::GetRandomEnergySR(G4double gamma,
                                                        G4double perpB)
{
  G4int i;
  static constexpr G4int iMax = 200;
  G4double energySR, random, position;

  random = G4UniformRand();

  for(i = 0; i < iMax; ++i)
  {
    if(random >= fIntegralProbabilityOfSR[i])
      break;
  }
  if(i <= 0)
    position = G4UniformRand();
  else if(i >= iMax)
    position = G4double(iMax);
  else
    position = i + G4UniformRand();

  energySR =
    0.0001 * position * position * fEnergyConst * gamma * gamma * perpB;

  if(energySR < 0.)
    energySR = 0.;

  return energySR;
}

/////////////////////////////////////////////////////////////////////////
G4double G4SynchrotronRadiationInMat::GetProbSpectrumSRforInt(G4double t)
{
  G4double result, hypCos2, hypCos = std::cosh(t);

  hypCos2 = hypCos * hypCos;
  result = std::cosh(5. * t / 3.) * std::exp(t - fKsi * hypCos);  // fKsi > 0. !
  result /= hypCos2;
  return result;
}

///////////////////////////////////////////////////////////////////////////
// return the probability to emit SR photon with relative energy
// energy/energy_c >= ksi
// for ksi <= 0. P = 1., however the method works for ksi > 0 only!
G4double G4SynchrotronRadiationInMat::GetIntProbSR(G4double ksi)
{
  if(ksi <= 0.)
    return 1.0;
  fKsi = ksi;  // should be > 0. !
  G4int n;
  G4double result, a;

  a = fAlpha;       // always = 0.
  n = fRootNumber;  // around default = 80

  G4Integrator<G4SynchrotronRadiationInMat,
               G4double (G4SynchrotronRadiationInMat::*)(G4double)>
    integral;

  result = integral.Laguerre(
    this, &G4SynchrotronRadiationInMat::GetProbSpectrumSRforInt, a, n);

  result *= 3. / 5. / pi;

  return result;
}

/////////////////////////////////////////////////////////////////////////
// return an auxiliary function for K_5/3 integral representation
G4double G4SynchrotronRadiationInMat::GetProbSpectrumSRforEnergy(G4double t)
{
  G4double result, hypCos = std::cosh(t);

  result = std::cosh(5. * t / 3.) * std::exp(t - fKsi * hypCos);  // fKsi > 0. !
  result /= hypCos;
  return result;
}

///////////////////////////////////////////////////////////////////////////
// return the probability to emit SR photon energy with relative energy
// energy/energy_c >= ksi
// for ksi <= 0. P = 1., however the method works for ksi > 0 only!
G4double G4SynchrotronRadiationInMat::GetEnergyProbSR(G4double ksi)
{
  if(ksi <= 0.)
    return 1.0;
  fKsi = ksi;  // should be > 0. !
  G4int n;
  G4double result, a;

  a = fAlpha;       // always = 0.
  n = fRootNumber;  // around default = 80

  G4Integrator<G4SynchrotronRadiationInMat,
               G4double (G4SynchrotronRadiationInMat::*)(G4double)>
    integral;

  result = integral.Laguerre(
    this, &G4SynchrotronRadiationInMat::GetProbSpectrumSRforEnergy, a, n);

  result *= 9. * std::sqrt(3.) * ksi / 8. / pi;

  return result;
}

/////////////////////////////////////////////////////////////////////////////
G4double G4SynchrotronRadiationInMat::GetIntegrandForAngleK(G4double t)
{
  G4double result, hypCos = std::cosh(t);

  result =
    std::cosh(fOrderAngleK * t) * std::exp(t - fEta * hypCos);  // fEta > 0. !
  result /= hypCos;
  return result;
}

//////////////////////////////////////////////////////////////////////////
// Return K 1/3 or 2/3 for angular distribution
G4double G4SynchrotronRadiationInMat::GetAngleK(G4double eta)
{
  fEta = eta;  // should be > 0. !
  G4int n;
  G4double result, a;

  a = fAlpha;       // always = 0.
  n = fRootNumber;  // around default = 80

  G4Integrator<G4SynchrotronRadiationInMat,
               G4double (G4SynchrotronRadiationInMat::*)(G4double)>
    integral;

  result = integral.Laguerre(
    this, &G4SynchrotronRadiationInMat::GetIntegrandForAngleK, a, n);

  return result;
}

/////////////////////////////////////////////////////////////////////////
// Relative angle diff distribution for given fKsi, which is set externally
G4double G4SynchrotronRadiationInMat::GetAngleNumberAtGammaKsi(G4double gpsi)
{
  G4double result, funK, funK2, gpsi2 = gpsi * gpsi;

  fPsiGamma = gpsi;
  fEta      = 0.5 * fKsi * (1. + gpsi2) * std::sqrt(1. + gpsi2);

  fOrderAngleK = 1. / 3.;
  funK         = GetAngleK(fEta);
  funK2        = funK * funK;

  result = gpsi2 * funK2 / (1. + gpsi2);

  fOrderAngleK = 2. / 3.;
  funK         = GetAngleK(fEta);
  funK2        = funK * funK;

  result += funK2;
  result *= (1. + gpsi2) * fKsi;

  return result;
}
