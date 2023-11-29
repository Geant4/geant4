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
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpRayleigh.cc
// Description: Discrete Process -- Rayleigh scattering of optical
//		photons
// Version:     1.0
// Created:     1996-05-31
// Author:      Juliet Armstrong
// Updated:     2014-10-10 -  This version calculates the Rayleigh scattering
//              length for more materials than just Water (although the Water
//              default is kept). To do this the user would need to specify the
//              ISOTHERMAL_COMPRESSIBILITY as a material property and
//              optionally an RS_SCALE_LENGTH (useful for testing). Code comes
//              from Philip Graham (Queen Mary University of London).
//              2010-06-11 - Fix Bug 207; Thanks to Xin Qian
//              (Kellogg Radiation Lab of Caltech)
//              2005-07-28 - add G4ProcessType to constructor
//              2001-10-18 by Peter Gumplinger
//              eliminate unused variable warning on Linux (gcc-2.95.2)
//              2001-09-18 by mma
//	          	>numOfMaterials=G4Material::GetNumberOfMaterials() in BuildPhy
//              2001-01-30 by Peter Gumplinger
//              > allow for positiv and negative CosTheta and force the
//              > new momentum direction to be in the same plane as the
//              > new and old polarization vectors
//              2001-01-29 by Peter Gumplinger
//              > fix calculation of SinTheta (from CosTheta)
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//
////////////////////////////////////////////////////////////////////////

#include "G4OpRayleigh.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpticalParameters.hh"
#include "G4OpProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpRayleigh::G4OpRayleigh(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  Initialise();
  SetProcessSubType(fOpRayleigh);
  thePhysicsTable = nullptr;

  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpRayleigh::~G4OpRayleigh()
{
  // VI: inside this PhysicsTable all properties are unique
  //     it is not possible to destroy
  delete thePhysicsTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpRayleigh::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpRayleigh::Initialise()
{
  SetVerboseLevel(G4OpticalParameters::Instance()->GetRayleighVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4OpRayleigh::PostStepDoIt(const G4Track& aTrack,
                                              const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

  if(verboseLevel > 1)
  {
    G4cout << "OpRayleigh: Scattering Photon!" << G4endl
           << "Old Momentum Direction: " << aParticle->GetMomentumDirection()
           << G4endl << "Old Polarization: " << aParticle->GetPolarization()
           << G4endl;
  }

  G4double cosTheta;
  G4ThreeVector oldMomDir, newMomDir;
  G4ThreeVector oldPol, newPol;
  G4double rand;
  G4double cost, sint, sinphi, cosphi;

  do
  {
    // Try to simulate the scattered photon momentum direction
    // w.r.t. the initial photon momentum direction
    cost = G4UniformRand();
    sint = std::sqrt(1. - cost * cost);
    // consider for the angle 90-180 degrees
    if(G4UniformRand() < 0.5)
      cost = -cost;

    // simulate the phi angle
    rand   = twopi * G4UniformRand();
    sinphi = std::sin(rand);
    cosphi = std::cos(rand);

    // construct the new momentum direction
    newMomDir.set(sint * cosphi, sint * sinphi, cost);
    oldMomDir = aParticle->GetMomentumDirection();
    newMomDir.rotateUz(oldMomDir);

    // calculate the new polarization direction
    // The new polarization needs to be in the same plane as the new
    // momentum direction and the old polarization direction
    oldPol = aParticle->GetPolarization();
    newPol = (oldPol - newMomDir.dot(oldPol) * newMomDir).unit();

    // There is a corner case, where the new momentum direction
    // is the same as old polarization direction:
    // random generate the azimuthal angle w.r.t. new momentum direction
    if(newPol.mag() == 0.)
    {
      rand = G4UniformRand() * twopi;
      newPol.set(std::cos(rand), std::sin(rand), 0.);
      newPol.rotateUz(newMomDir);
    }
    else
    {
      // There are two directions perpendicular to the new momentum direction
      if(G4UniformRand() < 0.5)
        newPol = -newPol;
    }

    // simulate according to the distribution cos^2(theta)
    cosTheta = newPol.dot(oldPol);
    // Loop checking, 13-Aug-2015, Peter Gumplinger
  } while(std::pow(cosTheta, 2) < G4UniformRand());

  aParticleChange.ProposePolarization(newPol);
  aParticleChange.ProposeMomentumDirection(newMomDir);

  if(verboseLevel > 1)
  {
    G4cout << "New Polarization: " << newPol << G4endl
           << "Polarization Change: " << *(aParticleChange.GetPolarization())
           << G4endl << "New Momentum Direction: " << newMomDir << G4endl
           << "Momentum Change: " << *(aParticleChange.GetMomentumDirection())
           << G4endl;
  }

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpRayleigh::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if(thePhysicsTable)
  {
    // thePhysicsTable->clearAndDestroy();
    delete thePhysicsTable;
    thePhysicsTable = nullptr;
  }

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const size_t numOfMaterials             = G4Material::GetNumberOfMaterials();
  thePhysicsTable                         = new G4PhysicsTable(numOfMaterials);

  for(size_t i = 0; i < numOfMaterials; ++i)
  {
    G4Material* material               = (*theMaterialTable)[i];
    G4MaterialPropertiesTable* matProp = material->GetMaterialPropertiesTable();
    G4PhysicsFreeVector* rayleigh = nullptr;
    if(matProp)
    {
      rayleigh = matProp->GetProperty(kRAYLEIGH);
      if(rayleigh == nullptr)
        rayleigh = CalculateRayleighMeanFreePaths(material);
    }
    thePhysicsTable->insertAt(i, rayleigh);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4OpRayleigh::GetMeanFreePath(const G4Track& aTrack, G4double,
                                       G4ForceCondition*)
{
  auto rayleigh = static_cast<G4PhysicsFreeVector*>(
      (*thePhysicsTable)(aTrack.GetMaterial()->GetIndex()));

  G4double rsLength = DBL_MAX;
  if(rayleigh)
  {
    rsLength = rayleigh->Value(aTrack.GetDynamicParticle()->GetTotalMomentum(),
                               idx_rslength);
  }
  return rsLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PhysicsFreeVector* G4OpRayleigh::CalculateRayleighMeanFreePaths(
  const G4Material* material) const
{
  G4MaterialPropertiesTable* MPT = material->GetMaterialPropertiesTable();

  // Retrieve the beta_T or isothermal compressibility value. For backwards
  // compatibility use a constant if the material is "Water". If the material
  // doesn't have an ISOTHERMAL_COMPRESSIBILITY constant then return
  G4double betat;
  if(material->GetName() == "Water")
  {
    betat = 7.658e-23 * m3 / MeV;
  }
  else if(MPT->ConstPropertyExists(kISOTHERMAL_COMPRESSIBILITY))
  {
    betat = MPT->GetConstProperty(kISOTHERMAL_COMPRESSIBILITY);
  }
  else
  {
    return nullptr;
  }

  // If the material doesn't have a RINDEX property vector then return
  G4MaterialPropertyVector* rIndex = MPT->GetProperty(kRINDEX);
  if(rIndex == nullptr)
    return nullptr;

  // Retrieve the optional scale factor (scales the scattering length)
  G4double scaleFactor = 1.0;
  if(MPT->ConstPropertyExists(kRS_SCALE_FACTOR))
  {
    scaleFactor = MPT->GetConstProperty(kRS_SCALE_FACTOR);
  }

  // Retrieve the material temperature. For backwards compatibility use a
  // constant if the material is "Water"
  G4double temperature;
  if(material->GetName() == "Water")
  {
    temperature =
      283.15 * kelvin;  // Temperature of water is 10 degrees celsius
  }
  else
  {
    temperature = material->GetTemperature();
  }

  auto rayleighMFPs = new G4PhysicsFreeVector();
  // This calculates the meanFreePath via the Einstein-Smoluchowski formula
  const G4double c1 =
    scaleFactor * betat * temperature * k_Boltzmann / (6.0 * pi);

  for(size_t uRIndex = 0; uRIndex < rIndex->GetVectorLength(); ++uRIndex)
  {
    const G4double energy        = rIndex->Energy(uRIndex);
    const G4double rIndexSquared = (*rIndex)[uRIndex] * (*rIndex)[uRIndex];
    const G4double xlambda       = h_Planck * c_light / energy;
    const G4double c2            = std::pow(twopi / xlambda, 4);
    const G4double c3 =
      std::pow(((rIndexSquared - 1.0) * (rIndexSquared + 2.0) / 3.0), 2);

    const G4double meanFreePath = 1.0 / (c1 * c2 * c3);

    if(verboseLevel > 0)
    {
      G4cout << energy << "MeV\t" << meanFreePath << "mm" << G4endl;
    }

    rayleighMFPs->InsertValues(energy, meanFreePath);
  }

  return rayleighMFPs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpRayleigh::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetRayleighVerboseLevel(verboseLevel);
}
