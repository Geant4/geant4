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
#include "G4QuasiCerenkov.hh"
#include "G4CerenkovQuasiTrackInfo.hh"

#include "G4ios.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleMomentum.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4Poisson.hh"
#include "G4QuasiOpticalData.hh"
#include "G4QuasiOpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4QuasiCerenkov::G4QuasiCerenkov(const G4String& processName, G4ProcessType type)
  : G4VProcess(processName, type)
  , fNumPhotons(0)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_QuasiCerenkov");
  SetProcessSubType(fCerenkov);

  thePhysicsTable = nullptr;

  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created." << G4endl;
  }
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4QuasiCerenkov::~G4QuasiCerenkov()
{
  if(thePhysicsTable != nullptr)
  {
    thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
  }
}

void G4QuasiCerenkov::ProcessDescription(std::ostream& out) const
{
  out << "The Cerenkov effect simulates optical photons created by the\n";
  out << "passage of charged particles through matter. Materials need\n";
  out << "to have the property RINDEX (refractive index) defined.\n";
  G4VProcess::DumpInfo();

  G4OpticalParameters* params = G4OpticalParameters::Instance();
  out << "Maximum beta change per step: " << params->GetCerenkovMaxBetaChange();
  out << "Maximum photons per step: " << params->GetCerenkovMaxPhotonsPerStep();
  out << "Track secondaries first: "
      << params->GetCerenkovTrackSecondariesFirst();
  out << "Stack photons: " << params->GetCerenkovStackPhotons();
  out << "Verbose level: " << params->GetCerenkovVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4QuasiCerenkov::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return (aParticleType.GetPDGCharge() != 0.0 &&
          aParticleType.GetPDGMass() != 0.0 &&
          aParticleType.GetParticleName() != "chargedgeantino" &&
          !aParticleType.IsShortLived())
           ? true
           : false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::Initialise()
{
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetMaxBetaChangePerStep(params->GetCerenkovMaxBetaChange());
  SetMaxNumPhotonsPerStep(params->GetCerenkovMaxPhotonsPerStep());
  SetTrackSecondariesFirst(params->GetCerenkovTrackSecondariesFirst());
  SetStackPhotons(params->GetCerenkovStackPhotons());
  SetOffloadPhotons(params->GetCerenkovOffloadPhotons());
  SetVerboseLevel(params->GetCerenkovVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(thePhysicsTable)
    return;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  std::size_t numOfMaterials              = G4Material::GetNumberOfMaterials();

  // Find the number of materials that have non-empty material property tables
  std::size_t numOfMaterialsWithMPT = 0;
  for(std::size_t i = 0; i < numOfMaterials; ++i)
  {
    if(((*theMaterialTable)[i])->GetMaterialPropertiesTable())
    {
       ++numOfMaterialsWithMPT;
    }
  }

  thePhysicsTable = new G4PhysicsTable(numOfMaterialsWithMPT);

  // loop over materials
  std::size_t indexMPT = 0;
  for(std::size_t i = 0; i < numOfMaterials; ++i)
  {
    G4PhysicsFreeVector* cerenkovIntegral = nullptr;

    // Retrieve vector of refraction indices for the material
    // from the material's optical properties table
    G4Material* aMaterial          = (*theMaterialTable)[i];
    G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();

    if(MPT)
    {
      cerenkovIntegral                          = new G4PhysicsFreeVector();
      G4MaterialPropertyVector* refractiveIndex = MPT->GetProperty(kRINDEX);

      if(refractiveIndex)
      {
        // Retrieve the first refraction index in vector
        // of (photon energy, refraction index) pairs
        G4double currentRI = (*refractiveIndex)[0];
        if(currentRI > 1.0)
        {
          // Create first (photon energy, Cerenkov Integral) pair
          G4double currentPM  = refractiveIndex->Energy(0);
          G4double currentCAI = 0.0;

          cerenkovIntegral->InsertValues(currentPM, currentCAI);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCAI = currentCAI;
          G4double prevRI  = currentRI;

          // loop over all (photon energy, refraction index)
          // pairs stored for this material
          for(std::size_t ii = 1; ii < refractiveIndex->GetVectorLength(); ++ii)
          {
            currentRI  = (*refractiveIndex)[ii];
            currentPM  = refractiveIndex->Energy(ii);
            currentCAI = prevCAI + (currentPM - prevPM) * 0.5 *
                                     (1.0 / (prevRI * prevRI) +
                                      1.0 / (currentRI * currentRI));

            cerenkovIntegral->InsertValues(currentPM, currentCAI);

            prevPM  = currentPM;
            prevCAI = currentCAI;
            prevRI  = currentRI;
          }
        }
      }
      // The Cerenkov integral for a given material will be inserted in
      // thePhysicsTable according to the position of the material in
      // the material table.
      thePhysicsTable->insertAt(indexMPT, cerenkovIntegral);
      fIndexMPT.insert(std::make_pair(i, indexMPT));
      ++indexMPT;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4QuasiCerenkov::PostStepDoIt(const G4Track& aTrack,
                                            const G4Step& aStep)
// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The
// parameters are then transformed into the Master Reference System, and
// they are added to the particle change.

{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial        = aTrack.GetMaterial();

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0      = pPreStepPoint->GetGlobalTime();

  G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
  if(!MPT)
    return pParticleChange;

  G4MaterialPropertyVector* Rindex = MPT->GetProperty(kRINDEX);
  if(!Rindex)
    return pParticleChange;

  G4double charge = aParticle->GetDefinition()->GetPDGCharge();

  G4double beta1 = pPreStepPoint->GetBeta();
  G4double beta2 = pPostStepPoint->GetBeta();
  G4double beta = (beta1 + beta2) * 0.5;

  G4double MeanNumberOfPhotons =
    GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);
  G4double MeanNumberOfPhotons1 =
    GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
  G4double MeanNumberOfPhotons2 =
    GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);

  if(MeanNumberOfPhotons <= 0.0)
  {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  MeanNumberOfPhotons *= aStep.GetStepLength();
  fNumPhotons         = (G4int) G4Poisson(MeanNumberOfPhotons);

  // third condition added to prevent infinite loop in do-while below,
  // see bugzilla 2555
  if(fNumPhotons <= 0 || !fStackingFlag ||
     std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2) < 1e-15)
  {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return pParticleChange;
  }

  G4double deltaVelocity = pPostStepPoint->GetVelocity() -
    pPreStepPoint->GetVelocity();
  auto touchableHandle = aStep.GetPreStepPoint()->GetTouchableHandle();

  if(fOffloadingFlag)
  {
    // Create a G4DynamicParticle with G4QuasiOpticalPhoton
    auto quasiPhoton = new G4DynamicParticle(
      G4QuasiOpticalPhoton::QuasiOpticalPhotonDefinition(),
      aParticle->GetMomentum());

    // Create a new G4Track object with the quasi-optical photon
    G4Track* aSecondaryTrack = new G4Track(quasiPhoton, t0, x0);
    aSecondaryTrack->SetTouchableHandle(touchableHandle);
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    aSecondaryTrack->SetCreatorModelID(secID);

    // Attach auxiliary track information with associated metadata
    G4QuasiOpticalData quasiTrackData{aMaterial->GetIndex(), fNumPhotons,
      charge, aStep.GetStepLength(), pPreStepPoint->GetVelocity(),
      deltaVelocity, aStep.GetDeltaPosition()};
    aSecondaryTrack->SetAuxiliaryTrackInformation(secID,
      new G4CerenkovQuasiTrackInfo(quasiTrackData,
        MeanNumberOfPhotons1, MeanNumberOfPhotons2));
    aParticleChange.AddSecondary(aSecondaryTrack);

    // Return early when offloading is enabled
    return pParticleChange;
  }

  ////////////////////////////////////////////////////////////////
  aParticleChange.SetNumberOfSecondaries(fNumPhotons);

  if(fTrackSecondariesFirst)
  {
    if(aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  }

  ////////////////////////////////////////////////////////////////
  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();
  G4double dp   = Pmax - Pmin;
  G4double deltaNumberOfPhotons = MeanNumberOfPhotons1 - MeanNumberOfPhotons2;
  G4double maxNumberOfPhotons =
    std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2);

  G4double nMax        = Rindex->GetMaxValue();
  G4double BetaInverse = 1. / beta;

  G4double maxCos  = BetaInverse / nMax;
  G4double maxSin2 = 1.0 - maxCos * maxCos;

  for(G4int i = 0; i < fNumPhotons; ++i)
  {
    // Determine photon energy
    G4double rand;
    G4double sampledEnergy;
    G4double cosTheta, sin2Theta;

    // sample an energy
    do
    {
      rand          = G4UniformRand();
      sampledEnergy = Pmin + rand * dp;
      cosTheta      = BetaInverse / Rindex->Value(sampledEnergy);

      sin2Theta = 1.0 - cosTheta * cosTheta;
      rand      = G4UniformRand();

      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while(rand * maxSin2 > sin2Theta);

    // Create photon momentum direction vector. The momentum direction is still
    // with respect to the coordinate system where the primary particle
    // direction is aligned with the z axis
    rand              = G4UniformRand();
    G4double phi      = twopi * rand;
    G4double sinPhi   = std::sin(phi);
    G4double cosPhi   = std::cos(phi);
    G4double sinTheta = std::sqrt(sin2Theta);
    G4ParticleMomentum photonMomentum(sinTheta * cosPhi, sinTheta * sinPhi,
                                      cosTheta);

    // Rotate momentum direction back to global reference system
    photonMomentum.rotateUz(p0);

    // Determine polarization of new photon
    G4ThreeVector photonPolarization(cosTheta * cosPhi, cosTheta * sinPhi,
                                     -sinTheta);

    // Rotate back to original coord system
    photonPolarization.rotateUz(p0);

    // Generate a new photon:
    auto aCerenkovPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);

    aCerenkovPhoton->SetPolarization(photonPolarization);
    aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

    G4double NumberOfPhotons;

    do
    {
      rand            = G4UniformRand();
      NumberOfPhotons = MeanNumberOfPhotons1 - rand * deltaNumberOfPhotons;
      // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
    } while(G4UniformRand() * maxNumberOfPhotons > NumberOfPhotons);

    G4double delta = rand * aStep.GetStepLength();
    G4double deltaTime =
      delta / (pPreStepPoint->GetVelocity() + rand * deltaVelocity * 0.5);

    G4double aSecondaryTime          = t0 + deltaTime;
    G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

    // Generate new G4Track object:
    G4Track* aSecondaryTrack =
      new G4Track(aCerenkovPhoton, aSecondaryTime, aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle(touchableHandle);
    aSecondaryTrack->SetParentID(aTrack.GetTrackID());
    aSecondaryTrack->SetCreatorModelID(secID);
    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if(verboseLevel > 1)
  {
    G4cout << "\n Exiting from G4QuasiCerenkov::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return pParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4QuasiCerenkov::GetMeanFreePath(const G4Track&, G4double,
                                          G4ForceCondition*)
{
  return 1.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4QuasiCerenkov::PostStepGetPhysicalInteractionLength(
  const G4Track& aTrack, G4double, G4ForceCondition* condition)
{
  *condition         = NotForced;
  G4double StepLimit = DBL_MAX;
  fNumPhotons        = 0;

  const G4Material* aMaterial = aTrack.GetMaterial();
  std::size_t materialIndex   = aMaterial->GetIndex();

  // If Physics Vector is not defined no Cerenkov photons
  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  G4MaterialPropertiesTable* MPT =
    ((*materialTable)[materialIndex])->GetMaterialPropertiesTable();
  //  if(!(*thePhysicsTable)[materialIndex])
  if(!MPT)
  {
    return StepLimit;
  }

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();

  G4double kineticEnergy                   = aParticle->GetKineticEnergy();
  const G4ParticleDefinition* particleType = aParticle->GetDefinition();
  G4double mass                            = particleType->GetPDGMass();

  G4double beta  = aParticle->GetTotalMomentum() / aParticle->GetTotalEnergy();
  G4double gamma = aParticle->GetTotalEnergy() / mass;

  G4MaterialPropertiesTable* aMaterialPropertiesTable =
    aMaterial->GetMaterialPropertiesTable();

  G4MaterialPropertyVector* Rindex = nullptr;

  if(aMaterialPropertiesTable)
    Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);

  G4double nMax;
  if(Rindex)
  {
    nMax = Rindex->GetMaxValue();
  }
  else
  {
    return StepLimit;
  }

  G4double BetaMin = 1. / nMax;
  if(BetaMin >= 1.)
    return StepLimit;

  G4double GammaMin = 1. / std::sqrt(1. - BetaMin * BetaMin);
  if(gamma < GammaMin)
    return StepLimit;

  G4double kinEmin = mass * (GammaMin - 1.);
  G4double RangeMin =
    G4LossTableManager::Instance()->GetRange(particleType, kinEmin, couple);
  G4double Range = G4LossTableManager::Instance()->GetRange(
    particleType, kineticEnergy, couple);
  G4double Step = Range - RangeMin;

  // If the step is smaller than G4ThreeVector::getTolerance(), it may happen
  // that the particle does not move. See bug 1992.
  static const G4double minAllowedStep = G4ThreeVector::getTolerance();
  if(Step < minAllowedStep)
    return StepLimit;

  if(Step < StepLimit)
    StepLimit = Step;

  // If user has defined an average maximum number of photons to be generated in
  // a Step, then calculate the Step length for that number of photons.
  if(fMaxPhotons > 0)
  {
    const G4double charge = aParticle->GetDefinition()->GetPDGCharge();
    G4double MeanNumberOfPhotons =
      GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);
    Step = 0.;
    if(MeanNumberOfPhotons > 0.0)
      Step = fMaxPhotons / MeanNumberOfPhotons;
    if(Step > 0. && Step < StepLimit)
      StepLimit = Step;
  }

  // If user has defined an maximum allowed change in beta per step
  if(fMaxBetaChange > 0.)
  {
    G4double dedx = G4LossTableManager::Instance()->GetDEDX(
      particleType, kineticEnergy, couple);
    G4double deltaGamma =
      gamma - 1. / std::sqrt(1. - beta * beta * (1. - fMaxBetaChange) *
                                    (1. - fMaxBetaChange));

    Step = mass * deltaGamma / dedx;
    if(Step > 0. && Step < StepLimit)
      StepLimit = Step;
  }

  *condition = StronglyForced;
  return StepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4QuasiCerenkov::GetAverageNumberOfPhotons(
  const G4double charge, const G4double beta, const G4Material* aMaterial,
  G4MaterialPropertyVector* Rindex) const
// This routine computes the number of Cerenkov photons produced per
// Geant4-unit (millimeter) in the current medium.
{
  constexpr G4double Rfact = 369.81 / (eV * cm);
  if(beta <= 0.0)
    return 0.0;
  G4double BetaInverse = 1. / beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  // 	- Refraction Indices for the current material
  //	- new G4PhysicsFreeVector allocated to hold CAI's
  std::size_t materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material
  auto it = fIndexMPT.find(materialIndex);

  std::size_t indexMPT = 0;
  if(it != fIndexMPT.end())
  {
    indexMPT = it->second;
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "G4MaterialPropertiesTable for " << aMaterial->GetName()
       << " is not found!" << G4endl;
    G4Exception("G4QuasiCerenkov::GetAverageNumberOfPhotons",
		"QuasiCerenkovCerenkov01", FatalException, ed);
  }

  G4PhysicsVector* CerenkovAngleIntegrals = ((*thePhysicsTable)(indexMPT));

  std::size_t length = CerenkovAngleIntegrals->GetVectorLength();
  if(0 == length)
    return 0.0;

  // Min and Max photon energies
  G4double Pmin = Rindex->Energy(0);
  G4double Pmax = Rindex->GetMaxEnergy();

  // Min and Max Refraction Indices
  G4double nMin = Rindex->GetMinValue();
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral
  G4double CAImax = (*CerenkovAngleIntegrals)[length - 1];

  G4double dp, ge;
  // If n(Pmax) < 1/Beta -- no photons generated
  if(nMax < BetaInverse)
  {
    dp = 0.0;
    ge = 0.0;
  }
  // otherwise if n(Pmin) >= 1/Beta -- photons generated
  else if(nMin > BetaInverse)
  {
    dp = Pmax - Pmin;
    ge = CAImax;
  }
  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then we need to find a P such
  // that the value of n(P) == 1/Beta. Interpolation is performed by the
  // GetEnergy() and Value() methods of the G4MaterialPropertiesTable and
  // the Value() method of G4PhysicsVector.
  else
  {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp   = Pmax - Pmin;

    G4double CAImin = CerenkovAngleIntegrals->Value(Pmin);
    ge              = CAImax - CAImin;

    if(verboseLevel > 1)
    {
      G4cout << "CAImin = " << CAImin << G4endl << "ge = " << ge << G4endl;
    }
  }

  // Calculate number of photons
  G4double NumPhotons = Rfact * charge / eplus * charge / eplus *
                        (dp - ge * BetaInverse * BetaInverse);

  return NumPhotons;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(
    fTrackSecondariesFirst);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::SetMaxBetaChangePerStep(const G4double value)
{
  fMaxBetaChange = value * CLHEP::perCent;
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons)
{
  fMaxPhotons = NumPhotons;
  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(fMaxPhotons);
}

void G4QuasiCerenkov::SetStackPhotons(const G4bool stackingFlag)
{
  fStackingFlag = stackingFlag;
  G4OpticalParameters::Instance()->SetCerenkovStackPhotons(fStackingFlag);
}

void G4QuasiCerenkov::SetOffloadPhotons(const G4bool offloadingFlag)
{
  fOffloadingFlag = offloadingFlag;
  G4OpticalParameters::Instance()->SetCerenkovOffloadPhotons(fOffloadingFlag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::DumpPhysicsTable() const
{
  G4cout << "Dump Physics Table!" << G4endl;
  for(std::size_t i = 0; i < thePhysicsTable->entries(); ++i)
  {
    (*thePhysicsTable)[i]->DumpValues();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4QuasiCerenkov::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetCerenkovVerboseLevel(verboseLevel);
}
