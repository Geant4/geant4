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
////////////////////////////////////////////////////////////////////////
//
// File G4OpMieHG.hh
// Description: Discrete Process -- Mie Scattering of Optical Photons
// Created: 2010-07-03
// Author: Xin Qian
// Based on work from Vlasios Vasileiou
//
// This subroutine will mimic the Mie scattering based on
// Henyey-Greenstein phase function
// Forward and backward angles are treated separately.
//
////////////////////////////////////////////////////////////////////////

#include "G4OpMieHG.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpticalParameters.hh"
#include "G4OpProcessSubType.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpMieHG::G4OpMieHG(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  Initialise();
  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  SetProcessSubType(fOpMieHG);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpMieHG::~G4OpMieHG() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpMieHG::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpMieHG::Initialise()
{
  SetVerboseLevel(G4OpticalParameters::Instance()->GetMieVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4OpMieHG::PostStepDoIt(const G4Track& aTrack,
                                           const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  const G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();

  G4double forwardRatio = MPT->GetConstProperty(kMIEHG_FORWARD_RATIO);

  if(verboseLevel > 1)
  {
    G4cout << "OpMie Scattering Photon!" << G4endl
           << " Old Momentum Direction: " << aParticle->GetMomentumDirection()
           << G4endl
           << " MIE Old Polarization: " << aParticle->GetPolarization()
           << G4endl;
  }

  G4bool isForward = (G4UniformRand() <= forwardRatio);
  G4double gg =
    MPT->GetConstProperty(isForward ? kMIEHG_FORWARD : kMIEHG_BACKWARD);

  // Sample the direction
  G4double r = G4UniformRand();
  G4double costh;

  if(gg != 0.)
  {
    G4double fact = (1. + gg)/(1. - gg + 2. * gg * r);
    costh = 2. * r * fact * fact * (1. - gg + gg * r) - 1.;
  }
  else
  {
    costh = 2. * r - 1.;
  }
  G4double sinth = std::sqrt(std::max(0., (1. - costh) * ( 1. + costh)));

  if(!isForward)
    costh *= -1.0;  // backward scattering

  G4double phi = G4UniformRand() * twopi;
  G4ThreeVector newMomDir(sinth * std::cos(phi), sinth * std::sin(phi), costh);
  newMomDir.rotateUz(aParticle->GetMomentumDirection());
  newMomDir = newMomDir.unit();

  G4ThreeVector oldPol = aParticle->GetPolarization();
  G4ThreeVector newPol = (newMomDir - oldPol / newMomDir.dot(oldPol)).unit();

  if(newPol.mag() == 0.)
  {
    r = G4UniformRand() * twopi;
    newPol.set(std::cos(r), std::sin(r), 0.);
    newPol.rotateUz(newMomDir);
  }
  else
  {
    // There are two directions perpendicular to new momentum direction
    if(G4UniformRand() < 0.5)
      newPol = -newPol;
  }

  aParticleChange.ProposePolarization(newPol);
  aParticleChange.ProposeMomentumDirection(newMomDir);

  if(verboseLevel > 1)
  {
    G4cout << "OpMie New Polarization: " << newPol << G4endl
           << " Polarization Change: " << *(aParticleChange.GetPolarization())
           << G4endl << " New Momentum Direction: " << newMomDir << G4endl
           << " Momentum Change: " << *(aParticleChange.GetMomentumDirection())
           << G4endl;
  }

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4OpMieHG::GetMeanFreePath(const G4Track& aTrack, G4double,
                                    G4ForceCondition*)
{
  G4double attLength = DBL_MAX;
  G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();
  if(MPT)
  {
    G4MaterialPropertyVector* attVector = MPT->GetProperty(kMIEHG);
    if(attVector)
    {
      attLength = attVector->Value(
        aTrack.GetDynamicParticle()->GetTotalEnergy(), idx_mie);
    }
  }
  return attLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpMieHG::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetMieVerboseLevel(verboseLevel);
}
