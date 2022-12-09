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
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4Scintillation.cc
// Description: RestDiscrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     2010-10-20 Allow the scintillation yield to be a function
//              of energy deposited by particle type
//              Thanks to Zach Hartwig (Department of Nuclear
//              Science and Engineeering - MIT)
//              2010-09-22 by Peter Gumplinger
//              > scintillation rise time included, thanks to
//              > Martin Goettlich/DESY
//              2005-08-17 by Peter Gumplinger
//              > change variable name MeanNumPhotons -> MeanNumberOfPhotons
//              2005-07-28 by Peter Gumplinger
//              > add G4ProcessType to constructor
//              2004-08-05 by Peter Gumplinger
//              > changed StronglyForced back to Forced in GetMeanLifeTime
//              2002-11-21 by Peter Gumplinger
//              > change to use G4Poisson for small MeanNumberOfPhotons
//              2002-11-07 by Peter Gumplinger
//              > now allow for fast and slow scintillation component
//              2002-11-05 by Peter Gumplinger
//              > now use scintillation constants from G4Material
//              2002-05-09 by Peter Gumplinger
//              > use only the PostStepPoint location for the origin of
//                scintillation photons when energy is lost to the medium
//                by a neutral particle
//              2000-09-18 by Peter Gumplinger
//              > change: aSecondaryPosition=x0+rand*aStep.GetDeltaPosition();
//                        aSecondaryTrack->SetTouchable(0);
//              2001-09-17, migration of Materials to pure STL (mma)
//              2003-06-03, V.Ivanchenko fix compilation warnings
//
////////////////////////////////////////////////////////////////////////

#include "G4Scintillation.hh"

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4EmProcessSubType.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4OpticalParameters.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsTable.hh"
#include "G4Poisson.hh"
#include "G4ScintillationTrackInformation.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Scintillation::G4Scintillation(const G4String& processName,
                                 G4ProcessType type)
  : G4VRestDiscreteProcess(processName, type)
  , fIntegralTable1(nullptr)
  , fIntegralTable2(nullptr)
  , fIntegralTable3(nullptr)
  , fEmSaturation(nullptr)
  , fNumPhotons(0)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_Scintillation");
  SetProcessSubType(fScintillation);

#ifdef G4DEBUG_SCINTILLATION
  ScintTrackEDep  = 0.;
  ScintTrackYield = 0.;
#endif

  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Scintillation::~G4Scintillation()
{
  if(fIntegralTable1 != nullptr)
  {
    fIntegralTable1->clearAndDestroy();
    delete fIntegralTable1;
  }
  if(fIntegralTable2 != nullptr)
  {
    fIntegralTable2->clearAndDestroy();
    delete fIntegralTable2;
  }
  if(fIntegralTable3 != nullptr)
  {
    fIntegralTable3->clearAndDestroy();
    delete fIntegralTable3;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::ProcessDescription(std::ostream& out) const
{
  out << "Scintillation simulates production of optical photons produced\n"
         "by a high energy particle traversing matter.\n"
         "Various material properties need to be defined.\n";
  G4VRestDiscreteProcess::DumpInfo();

  G4OpticalParameters* params = G4OpticalParameters::Instance();
  out << "Track secondaries first: " << params->GetScintTrackSecondariesFirst();
  out << "Finite rise time: " << params->GetScintFiniteRiseTime();
  out << "Scintillation by particle type: " << params->GetScintByParticleType();
  out << "Save track information: " << params->GetScintTrackInfo();
  out << "Stack photons: " << params->GetScintStackPhotons();
  out << "Verbose level: " << params->GetScintVerboseLevel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4Scintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  if(aParticleType.GetParticleName() == "opticalphoton")
    return false;
  if(aParticleType.IsShortLived())
    return false;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::Initialise()
{
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetTrackSecondariesFirst(params->GetScintTrackSecondariesFirst());
  SetFiniteRiseTime(params->GetScintFiniteRiseTime());
  SetScintillationByParticleType(params->GetScintByParticleType());
  SetScintillationTrackInfo(params->GetScintTrackInfo());
  SetStackPhotons(params->GetScintStackPhotons());
  SetVerboseLevel(params->GetScintVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(fIntegralTable1 != nullptr)
  {
    fIntegralTable1->clearAndDestroy();
    delete fIntegralTable1;
    fIntegralTable1 = nullptr;
  }
  if(fIntegralTable2 != nullptr)
  {
    fIntegralTable2->clearAndDestroy();
    delete fIntegralTable2;
    fIntegralTable2 = nullptr;
  }
  if(fIntegralTable3 != nullptr)
  {
    fIntegralTable3->clearAndDestroy();
    delete fIntegralTable3;
    fIntegralTable3 = nullptr;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  std::size_t numOfMaterials           = G4Material::GetNumberOfMaterials();

  // create new physics table
  if(!fIntegralTable1)
    fIntegralTable1 = new G4PhysicsTable(numOfMaterials);
  if(!fIntegralTable2)
    fIntegralTable2 = new G4PhysicsTable(numOfMaterials);
  if(!fIntegralTable3)
    fIntegralTable3 = new G4PhysicsTable(numOfMaterials);

  for(std::size_t i = 0; i < numOfMaterials; ++i)
  {
    auto vector1 = new G4PhysicsFreeVector();
    auto vector2 = new G4PhysicsFreeVector();
    auto vector3 = new G4PhysicsFreeVector();

    // Retrieve vector of scintillation wavelength intensity for
    // the material from the material's optical properties table.
    G4MaterialPropertiesTable* MPT =
      ((*materialTable)[i])->GetMaterialPropertiesTable();

    if(MPT)
    {
      G4MaterialPropertyVector* MPV =
        MPT->GetProperty(kSCINTILLATIONCOMPONENT1);
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector1->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(std::size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);

            vector1->InsertValues(currentPM, currentCII);

            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }

      MPV = MPT->GetProperty(kSCINTILLATIONCOMPONENT2);
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector2->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(std::size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);

            vector2->InsertValues(currentPM, currentCII);

            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
      MPV = MPT->GetProperty(kSCINTILLATIONCOMPONENT3);
      if(MPV)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*MPV)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy, Scintillation Integral pair
          G4double currentPM  = MPV->Energy(0);
          G4double currentCII = 0.0;
          vector3->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(std::size_t ii = 1; ii < MPV->GetVectorLength(); ++ii)
          {
            currentPM = MPV->Energy(ii);
            currentIN = (*MPV)[ii];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);

            vector3->InsertValues(currentPM, currentCII);

            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
    }
    fIntegralTable1->insertAt(i, vector1);
    fIntegralTable2->insertAt(i, vector2);
    fIntegralTable3->insertAt(i, vector3);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4Scintillation::AtRestDoIt(const G4Track& aTrack,
                                               const G4Step& aStep)
// This routine simply calls the equivalent PostStepDoIt since all the
// necessary information resides in aStep.GetTotalEnergyDeposit()
{
  return G4Scintillation::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4Scintillation::PostStepDoIt(const G4Track& aTrack,
                                                 const G4Step& aStep)
// This routine is called for each tracking step of a charged particle
// in a scintillator. A Poisson/Gauss-distributed number of photons is
// generated according to the scintillation yield formula, distributed
// evenly along the track segment and uniformly into 4pi.
{
  aParticleChange.Initialize(aTrack);
  fNumPhotons = 0;

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4Material* aMaterial        = aTrack.GetMaterial();

  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0      = pPreStepPoint->GetGlobalTime();

  G4double TotalEnergyDeposit = aStep.GetTotalEnergyDeposit();

  G4MaterialPropertiesTable* MPT = aMaterial->GetMaterialPropertiesTable();
  if(!MPT)
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);

  G4int N_timeconstants = 1;

  if(MPT->GetProperty(kSCINTILLATIONCOMPONENT3))
    N_timeconstants = 3;
  else if(MPT->GetProperty(kSCINTILLATIONCOMPONENT2))
    N_timeconstants = 2;
  else if(!(MPT->GetProperty(kSCINTILLATIONCOMPONENT1)))
  {
    // no components were specified
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4double ResolutionScale = MPT->GetConstProperty(kRESOLUTIONSCALE);
  G4double MeanNumberOfPhotons;

  G4double yield1     = 0.;
  G4double yield2     = 0.;
  G4double yield3     = 0.;
  G4double sum_yields = 0.;

  if(fScintillationByParticleType)
  {
    MeanNumberOfPhotons = GetScintillationYieldByParticleType(
      aTrack, aStep, yield1, yield2, yield3);
  }
  else
  {
    yield1 = MPT->ConstPropertyExists(kSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kSCINTILLATIONYIELD3)
               : 0.;
    // The default linear scintillation process
    // Units: [# scintillation photons / MeV]
    MeanNumberOfPhotons = MPT->GetConstProperty(kSCINTILLATIONYIELD);
    // Birk's correction via fEmSaturation and specifying scintillation by
    // by particle type are physically mutually exclusive
    if(fEmSaturation)
      MeanNumberOfPhotons *=
        (fEmSaturation->VisibleEnergyDepositionAtAStep(&aStep));
    else
      MeanNumberOfPhotons *= TotalEnergyDeposit;
  }
  sum_yields = yield1 + yield2 + yield3;

  if(MeanNumberOfPhotons > 10.)
  {
    G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
    fNumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons, sigma) + 0.5);
  }
  else
  {
    fNumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
  }

  if(fNumPhotons <= 0 || !fStackingFlag)
  {
    // return unchanged particle and no secondaries
    aParticleChange.SetNumberOfSecondaries(0);
    return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  aParticleChange.SetNumberOfSecondaries(fNumPhotons);

  if(fTrackSecondariesFirst)
  {
    if(aTrack.GetTrackStatus() == fAlive)
      aParticleChange.ProposeTrackStatus(fSuspend);
  }

  G4int materialIndex = (G4int)aMaterial->GetIndex();

  // Retrieve the Scintillation Integral for this material
  // new G4PhysicsFreeVector allocated to hold CII's
  std::size_t numPhot                = fNumPhotons;
  G4double scintTime                 = 0.;
  G4double riseTime                  = 0.;
  G4PhysicsFreeVector* scintIntegral = nullptr;
  G4ScintillationType scintType      = Slow;

  for(G4int scnt = 0; scnt < N_timeconstants; ++scnt)
  {
    // if there is 1 time constant it is #1, etc.
    if(scnt == 0)
    {
      if(N_timeconstants == 1)
      {
        numPhot = fNumPhotons;
      }
      else
      {
        numPhot = yield1 / sum_yields * fNumPhotons;
      }
      scintTime = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty(kSCINTILLATIONRISETIME1);
      }
      scintType = Fast;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable1)(materialIndex));
    }
    else if(scnt == 1)
    {
      // to be consistent with old version (due to double->int conversion)
      if(N_timeconstants == 2)
      {
        numPhot = fNumPhotons - numPhot;
      }
      else
      {
        numPhot = yield2 / sum_yields * fNumPhotons;
      }
      scintTime = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty(kSCINTILLATIONRISETIME2);
      }
      scintType = Medium;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable2)(materialIndex));
    }
    else if(scnt == 2)
    {
      numPhot   = yield3 / sum_yields * fNumPhotons;
      scintTime = MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT3);
      if(fFiniteRiseTime)
      {
        riseTime = MPT->GetConstProperty(kSCINTILLATIONRISETIME3);
      }
      scintType = Slow;
      scintIntegral =
        (G4PhysicsFreeVector*) ((*fIntegralTable3)(materialIndex));
    }

    if(!scintIntegral)
      continue;

    G4double CIImax = scintIntegral->GetMaxValue();
    for(std::size_t i = 0; i < numPhot; ++i)
    {
      // Determine photon energy
      G4double CIIvalue      = G4UniformRand() * CIImax;
      G4double sampledEnergy = scintIntegral->GetEnergy(CIIvalue);

      if(verboseLevel > 1)
      {
        G4cout << "sampledEnergy = " << sampledEnergy << G4endl;
        G4cout << "CIIvalue =        " << CIIvalue << G4endl;
      }

      // Generate random photon direction
      G4double cost = 1. - 2. * G4UniformRand();
      G4double sint = std::sqrt((1. - cost) * (1. + cost));
      G4double phi  = twopi * G4UniformRand();
      G4double sinp = std::sin(phi);
      G4double cosp = std::cos(phi);
      G4ParticleMomentum photonMomentum(sint * cosp, sint * sinp, cost);

      // Determine polarization of new photon
      G4ThreeVector photonPolarization(cost * cosp, cost * sinp, -sint);
      G4ThreeVector perp = photonMomentum.cross(photonPolarization);
      phi                = twopi * G4UniformRand();
      sinp               = std::sin(phi);
      cosp               = std::cos(phi);
      photonPolarization = (cosp * photonPolarization + sinp * perp).unit();

      // Generate a new photon:
      auto scintPhoton = new G4DynamicParticle(opticalphoton, photonMomentum);
      scintPhoton->SetPolarization(photonPolarization);
      scintPhoton->SetKineticEnergy(sampledEnergy);

      // Generate new G4Track object:
      G4double rand = G4UniformRand();
      if(aParticle->GetDefinition()->GetPDGCharge() == 0)
      {
        rand = 1.0;
      }

      // emission time distribution
      G4double delta = rand * aStep.GetStepLength();
      G4double deltaTime =
        delta /
        (pPreStepPoint->GetVelocity() +
         rand * (pPostStepPoint->GetVelocity() - pPreStepPoint->GetVelocity()) /
           2.);
      if(riseTime == 0.0)
      {
        deltaTime -= scintTime * std::log(G4UniformRand());
      }
      else
      {
        deltaTime += sample_time(riseTime, scintTime);
      }

      G4double secTime          = t0 + deltaTime;
      G4ThreeVector secPosition = x0 + rand * aStep.GetDeltaPosition();

      G4Track* secTrack = new G4Track(scintPhoton, secTime, secPosition);
      secTrack->SetTouchableHandle(
        aStep.GetPreStepPoint()->GetTouchableHandle());
      secTrack->SetParentID(aTrack.GetTrackID());
      secTrack->SetCreatorModelID(secID);
      if(fScintillationTrackInfo)
        secTrack->SetUserInformation(
          new G4ScintillationTrackInformation(scintType));
      aParticleChange.AddSecondary(secTrack);
    }
  }

  if(verboseLevel > 1)
  {
    G4cout << "\n Exiting from G4Scintillation::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return G4VRestDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4Scintillation::GetMeanFreePath(const G4Track&, G4double,
                                          G4ForceCondition* condition)
{
  *condition = StronglyForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4Scintillation::GetMeanLifeTime(const G4Track&,
                                          G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4Scintillation::sample_time(G4double tau1, G4double tau2)
{
  // tau1: rise time and tau2: decay time
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while(true)
  {
    G4double ran1 = G4UniformRand();
    G4double ran2 = G4UniformRand();

    // exponential distribution as envelope function: very efficient
    G4double d = (tau1 + tau2) / tau2;
    // make sure the envelope function is
    // always larger than the bi-exponential
    G4double t  = -1.0 * tau2 * std::log(1. - ran1);
    G4double gg = d * single_exp(t, tau2);
    if(ran2 <= bi_exp(t, tau1, tau2) / gg)
      return t;
  }
  return -1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4Scintillation::GetScintillationYieldByParticleType(
  const G4Track& aTrack, const G4Step& aStep, G4double& yield1,
  G4double& yield2, G4double& yield3)
{
  // new in 10.7, allow multiple time constants with ScintByParticleType
  // Get the G4MaterialPropertyVector containing the scintillation
  // yield as a function of the energy deposited and particle type

  G4ParticleDefinition* pDef = aTrack.GetDynamicParticle()->GetDefinition();
  G4MaterialPropertyVector* yieldVector = nullptr;
  G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();

  // Protons
  if(pDef == G4Proton::ProtonDefinition())
  {
    yieldVector = MPT->GetProperty(kPROTONSCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kPROTONSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kPROTONSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kPROTONSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kPROTONSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kPROTONSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kPROTONSCINTILLATIONYIELD3)
               : 0.;
  }

  // Deuterons
  else if(pDef == G4Deuteron::DeuteronDefinition())
  {
    yieldVector = MPT->GetProperty(kDEUTERONSCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kDEUTERONSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kDEUTERONSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kDEUTERONSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kDEUTERONSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kDEUTERONSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kDEUTERONSCINTILLATIONYIELD3)
               : 0.;
  }

  // Tritons
  else if(pDef == G4Triton::TritonDefinition())
  {
    yieldVector = MPT->GetProperty(kTRITONSCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kTRITONSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kTRITONSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kTRITONSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kTRITONSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kTRITONSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kTRITONSCINTILLATIONYIELD3)
               : 0.;
  }

  // Alphas
  else if(pDef == G4Alpha::AlphaDefinition())
  {
    yieldVector = MPT->GetProperty(kALPHASCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kALPHASCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kALPHASCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kALPHASCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kALPHASCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kALPHASCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kALPHASCINTILLATIONYIELD3)
               : 0.;
  }

  // Ions (particles derived from G4VIon and G4Ions) and recoil ions
  // below the production cut from neutrons after hElastic
  else if(pDef->GetParticleType() == "nucleus" ||
          pDef == G4Neutron::NeutronDefinition())
  {
    yieldVector = MPT->GetProperty(kIONSCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kIONSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kIONSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kIONSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kIONSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kIONSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kIONSCINTILLATIONYIELD3)
               : 0.;
  }

  // Electrons (must also account for shell-binding energy
  // attributed to gamma from standard photoelectric effect)
  // and, default for particles not enumerated/listed above
  else
  {
    yieldVector = MPT->GetProperty(kELECTRONSCINTILLATIONYIELD);
    yield1      = MPT->ConstPropertyExists(kELECTRONSCINTILLATIONYIELD1)
               ? MPT->GetConstProperty(kELECTRONSCINTILLATIONYIELD1)
               : 1.;
    yield2 = MPT->ConstPropertyExists(kELECTRONSCINTILLATIONYIELD2)
               ? MPT->GetConstProperty(kELECTRONSCINTILLATIONYIELD2)
               : 0.;
    yield3 = MPT->ConstPropertyExists(kELECTRONSCINTILLATIONYIELD3)
               ? MPT->GetConstProperty(kELECTRONSCINTILLATIONYIELD3)
               : 0.;
  }

  // Throw an exception if no scintillation yield vector is found
  if(!yieldVector)
  {
    G4ExceptionDescription ed;
    ed << "\nG4Scintillation::PostStepDoIt(): "
       << "Request for scintillation yield for energy deposit and particle\n"
       << "type without correct entry in MaterialPropertiesTable.\n"
       << "ScintillationByParticleType requires at minimum that \n"
       << "ELECTRONSCINTILLATIONYIELD is set by the user\n"
       << G4endl;
    G4String comments = "Missing MaterialPropertiesTable entry - No correct "
                        "entry in MaterialPropertiesTable";
    G4Exception("G4Scintillation::PostStepDoIt", "Scint01", FatalException, ed,
                comments);
  }

  ///////////////////////////////////////
  // Calculate the scintillation light //
  ///////////////////////////////////////
  // To account for potential nonlinearity and scintillation photon
  // density along the track, light (L) is produced according to:
  // L_currentStep = L(PreStepKE) - L(PreStepKE - EDep)

  G4double ScintillationYield   = 0.;
  G4double StepEnergyDeposit    = aStep.GetTotalEnergyDeposit();
  G4double PreStepKineticEnergy = aStep.GetPreStepPoint()->GetKineticEnergy();

  if(PreStepKineticEnergy <= yieldVector->GetMaxEnergy())
  {
    // G4double Yield1 = yieldVector->Value(PreStepKineticEnergy);
    // G4double Yield2 = yieldVector->Value(PreStepKineticEnergy -
    // StepEnergyDeposit); ScintillationYield = Yield1 - Yield2;
    ScintillationYield =
      yieldVector->Value(PreStepKineticEnergy) -
      yieldVector->Value(PreStepKineticEnergy - StepEnergyDeposit);
  }
  else
  {
    G4ExceptionDescription ed;
    ed << "\nG4Scintillation::GetScintillationYieldByParticleType(): Request\n"
       << "for scintillation light yield above the available energy range\n"
       << "specified in G4MaterialPropertiesTable. A linear interpolation\n"
       << "will be performed to compute the scintillation light yield using\n"
       << "(L_max / E_max) as the photon yield per unit energy." << G4endl;
    G4String cmt = "\nScintillation yield may be unphysical!\n";
    G4Exception("G4Scintillation::GetScintillationYieldByParticleType()",
                "Scint03", JustWarning, ed, cmt);

    // Units: [# scintillation photons]
    ScintillationYield = yieldVector->GetMaxValue() /
                         yieldVector->GetMaxEnergy() * StepEnergyDeposit;
  }

#ifdef G4DEBUG_SCINTILLATION
  // Increment track aggregators
  ScintTrackYield += ScintillationYield;
  ScintTrackEDep += StepEnergyDeposit;

  G4cout << "\n--- G4Scintillation::GetScintillationYieldByParticleType() ---\n"
         << "--\n"
         << "--  Name         =  "
         << aTrack.GetParticleDefinition()->GetParticleName() << "\n"
         << "--  TrackID      =  " << aTrack.GetTrackID() << "\n"
         << "--  ParentID     =  " << aTrack.GetParentID() << "\n"
         << "--  Current KE   =  " << aTrack.GetKineticEnergy() / MeV
         << " MeV\n"
         << "--  Step EDep    =  " << aStep.GetTotalEnergyDeposit() / MeV
         << " MeV\n"
         << "--  Track EDep   =  " << ScintTrackEDep / MeV << " MeV\n"
         << "--  Vertex KE    =  " << aTrack.GetVertexKineticEnergy() / MeV
         << " MeV\n"
         << "--  Step yield   =  " << ScintillationYield << " photons\n"
         << "--  Track yield  =  " << ScintTrackYield << " photons\n"
         << G4endl;

  // The track has terminated within or has left the scintillator volume
  if((aTrack.GetTrackStatus() == fStopButAlive) or
     (aStep.GetPostStepPoint()->GetStepStatus() == fGeomBoundary))
  {
    // Reset aggregators for the next track
    ScintTrackEDep  = 0.;
    ScintTrackYield = 0.;
  }
#endif

  return ScintillationYield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::DumpPhysicsTable() const
{
  if(fIntegralTable1)
  {
    for(std::size_t i = 0; i < fIntegralTable1->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable1)[i])->DumpValues();
    }
  }
  if(fIntegralTable2)
  {
    for(std::size_t i = 0; i < fIntegralTable2->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable2)[i])->DumpValues();
    }
  }
  if(fIntegralTable3)
  {
    for(std::size_t i = 0; i < fIntegralTable3->entries(); ++i)
    {
      ((G4PhysicsFreeVector*) (*fIntegralTable3)[i])->DumpValues();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetTrackSecondariesFirst(const G4bool state)
{
  fTrackSecondariesFirst = state;
  G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(
    fTrackSecondariesFirst);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetFiniteRiseTime(const G4bool state)
{
  fFiniteRiseTime = state;
  G4OpticalParameters::Instance()->SetScintFiniteRiseTime(fFiniteRiseTime);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetScintillationByParticleType(const G4bool scintType)
{
  if(fEmSaturation && scintType)
  {
    G4Exception("G4Scintillation::SetScintillationByParticleType", "Scint02",
                JustWarning,
                "Redefinition: Birks Saturation is replaced by "
                "ScintillationByParticleType!");
    RemoveSaturation();
  }
  fScintillationByParticleType = scintType;
  G4OpticalParameters::Instance()->SetScintByParticleType(
    fScintillationByParticleType);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetScintillationTrackInfo(const G4bool trackType)
{
  fScintillationTrackInfo = trackType;
  G4OpticalParameters::Instance()->SetScintTrackInfo(fScintillationTrackInfo);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetStackPhotons(const G4bool stackingFlag)
{
  fStackingFlag = stackingFlag;
  G4OpticalParameters::Instance()->SetScintStackPhotons(fStackingFlag);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4Scintillation::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetScintVerboseLevel(verboseLevel);
}
