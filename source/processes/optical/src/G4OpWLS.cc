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
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS.cc
// Description: Discrete Process -- Wavelength Shifting of Optical Photons
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
//
////////////////////////////////////////////////////////////////////////

#include "G4OpWLS.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"
#include "G4Poisson.hh"
#include "G4OpticalParameters.hh"
#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "G4WLSTimeGeneratorProfileExponential.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpWLS::G4OpWLS(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  WLSTimeGeneratorProfile = nullptr;
  Initialise();
  SetProcessSubType(fOpWLS);
  theIntegralTable = nullptr;

  if(verboseLevel > 0)
    G4cout << GetProcessName() << " is created " << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpWLS::~G4OpWLS()
{
  if(theIntegralTable)
  {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
  }
  delete WLSTimeGeneratorProfile;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpWLS::PreparePhysicsTable(const G4ParticleDefinition&) { Initialise(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpWLS::Initialise()
{
  G4OpticalParameters* params = G4OpticalParameters::Instance();
  SetVerboseLevel(params->GetWLSVerboseLevel());
  UseTimeProfile(params->GetWLSTimeProfile());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4OpWLS::PostStepDoIt(const G4Track& aTrack,
                                         const G4Step& aStep)
{
  std::vector<G4Track*> proposedSecondaries;
  aParticleChange.Initialize(aTrack);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if(verboseLevel > 1)
  {
    G4cout << "\n** G4OpWLS: Photon absorbed! **" << G4endl;
  }

  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();
  if(!MPT)
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  if(!MPT->GetProperty(kWLSCOMPONENT))
  {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4int NumPhotons = 1;
  if(MPT->ConstPropertyExists(kWLSMEANNUMBERPHOTONS))
  {
    G4double MeanNumberOfPhotons = MPT->GetConstProperty(kWLSMEANNUMBERPHOTONS);
    NumPhotons                   = G4int(G4Poisson(MeanNumberOfPhotons));
    if(NumPhotons <= 0)
    {
      // return unchanged particle and no secondaries
      aParticleChange.SetNumberOfSecondaries(0);
      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }

  // Retrieve the WLS Integral for this material
  // new G4PhysicsFreeVector allocated to hold CII's
  G4double primaryEnergy = aTrack.GetDynamicParticle()->GetKineticEnergy();
  G4double WLSTime       = 0.;
  G4PhysicsFreeVector* WLSIntegral = nullptr;

  WLSTime     = MPT->GetConstProperty(kWLSTIMECONSTANT);
  WLSIntegral = (G4PhysicsFreeVector*) ((*theIntegralTable)(
    aTrack.GetMaterial()->GetIndex()));

  // Max WLS Integral
  G4double CIImax       = WLSIntegral->GetMaxValue();
  G4int NumberOfPhotons = NumPhotons;

  for(G4int i = 0; i < NumPhotons; ++i)
  {
    G4double sampledEnergy;
    // Make sure the energy of the secondary is less than that of the primary
    for(G4int j = 1; j <= 100; ++j)
    {
      // Determine photon energy
      G4double CIIvalue = G4UniformRand() * CIImax;
      sampledEnergy     = WLSIntegral->GetEnergy(CIIvalue);
      if(sampledEnergy <= primaryEnergy)
        break;
    }
    // If no such energy can be sampled, return one less secondary, or none
    if(sampledEnergy > primaryEnergy)
    {
      if(verboseLevel > 1)
      {
        G4cout << " *** G4OpWLS: One less WLS photon will be returned ***"
               << G4endl;
      }
      NumberOfPhotons--;
      if(NumberOfPhotons == 0)
      {
        if(verboseLevel > 1)
        {
          G4cout
            << " *** G4OpWLS: No WLS photon can be sampled for this primary ***"
            << G4endl;
        }
        // return unchanged particle and no secondaries
        aParticleChange.SetNumberOfSecondaries(0);
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
      continue;
    }
    else if(verboseLevel > 1)
    {
      G4cout << "G4OpWLS: Created photon with energy: " << sampledEnergy
             << G4endl;
    }

    // Generate random photon direction
    G4double cost = 1. - 2. * G4UniformRand();
    G4double sint = std::sqrt((1. - cost) * (1. + cost));
    G4double phi  = twopi * G4UniformRand();
    G4double sinp = std::sin(phi);
    G4double cosp = std::cos(phi);
    G4ParticleMomentum photonMomentum(sint * cosp, sint * sinp, cost);

    G4ThreeVector photonPolarization(cost * cosp, cost * sinp, -sint);
    G4ThreeVector perp = photonMomentum.cross(photonPolarization);

    phi                = twopi * G4UniformRand();
    sinp               = std::sin(phi);
    cosp               = std::cos(phi);
    photonPolarization = (cosp * photonPolarization + sinp * perp).unit();

    // Generate a new photon:
    auto sec_dp =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
    sec_dp->SetPolarization(photonPolarization);
    sec_dp->SetKineticEnergy(sampledEnergy);

    G4double secTime = pPostStepPoint->GetGlobalTime() +
                       WLSTimeGeneratorProfile->GenerateTime(WLSTime);
    G4ThreeVector secPos = pPostStepPoint->GetPosition();
    G4Track* secTrack    = new G4Track(sec_dp, secTime, secPos);

    secTrack->SetTouchableHandle(aTrack.GetTouchableHandle());
    secTrack->SetParentID(aTrack.GetTrackID());

    proposedSecondaries.push_back(secTrack);
  }

  aParticleChange.SetNumberOfSecondaries((G4int)proposedSecondaries.size());
  for(auto sec : proposedSecondaries)
  {
    aParticleChange.AddSecondary(sec);
  }
  if(verboseLevel > 1)
  {
    G4cout << "\n Exiting from G4OpWLS::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpWLS::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(theIntegralTable)
  {
    theIntegralTable->clearAndDestroy();
    delete theIntegralTable;
    theIntegralTable = nullptr;
  }

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();
  std::size_t numOfMaterials           = G4Material::GetNumberOfMaterials();
  theIntegralTable                     = new G4PhysicsTable(numOfMaterials);

  // loop for materials
  for(std::size_t i = 0; i < numOfMaterials; ++i)
  {
    auto physVector = new G4PhysicsFreeVector();

    // Retrieve vector of WLS wavelength intensity for
    // the material from the material's optical properties table.
    G4MaterialPropertiesTable* MPT =
      (*materialTable)[i]->GetMaterialPropertiesTable();
    if(MPT)
    {
      G4MaterialPropertyVector* wlsVector = MPT->GetProperty(kWLSCOMPONENT);
      if(wlsVector)
      {
        // Retrieve the first intensity point in vector
        // of (photon energy, intensity) pairs
        G4double currentIN = (*wlsVector)[0];
        if(currentIN >= 0.0)
        {
          // Create first (photon energy)
          G4double currentPM  = wlsVector->Energy(0);
          G4double currentCII = 0.0;
          physVector->InsertValues(currentPM, currentCII);

          // Set previous values to current ones prior to loop
          G4double prevPM  = currentPM;
          G4double prevCII = currentCII;
          G4double prevIN  = currentIN;

          // loop over all (photon energy, intensity)
          // pairs stored for this material
          for(std::size_t j = 1; j < wlsVector->GetVectorLength(); ++j)
          {
            currentPM = wlsVector->Energy(j);
            currentIN = (*wlsVector)[j];
            currentCII =
              prevCII + 0.5 * (currentPM - prevPM) * (prevIN + currentIN);

            physVector->InsertValues(currentPM, currentCII);

            prevPM  = currentPM;
            prevCII = currentCII;
            prevIN  = currentIN;
          }
        }
      }
    }
    theIntegralTable->insertAt(i, physVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4OpWLS::GetMeanFreePath(const G4Track& aTrack, G4double,
                                  G4ForceCondition*)
{
  G4double thePhotonEnergy = aTrack.GetDynamicParticle()->GetTotalEnergy();
  G4double attLength       = DBL_MAX;
  G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();

  if(MPT)
  {
    G4MaterialPropertyVector* attVector = MPT->GetProperty(kWLSABSLENGTH);
    if(attVector)
    {
      attLength = attVector->Value(thePhotonEnergy, idx_wls);
    }
  }
  return attLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpWLS::UseTimeProfile(const G4String name)
{
  if(WLSTimeGeneratorProfile)
  {
    delete WLSTimeGeneratorProfile;
    WLSTimeGeneratorProfile = nullptr;
  }
  if(name == "delta")
  {
    WLSTimeGeneratorProfile = new G4WLSTimeGeneratorProfileDelta("delta");
  }
  else if(name == "exponential")
  {
    WLSTimeGeneratorProfile =
      new G4WLSTimeGeneratorProfileExponential("exponential");
  }
  else
  {
    G4Exception("G4OpWLS::UseTimeProfile", "em0202", FatalException,
                "generator does not exist");
  }
  G4OpticalParameters::Instance()->SetWLSTimeProfile(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpWLS::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetWLSVerboseLevel(verboseLevel);
}
