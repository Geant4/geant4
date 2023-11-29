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
// Optical Photon Absorption Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpAbsorption.cc
// Description: Discrete Process -- Absorption of Optical Photons
// Version:     1.0
// Created:     1996-05-21
// Author:      Juliet Armstrong
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2000-09-18 by Peter Gumplinger
//              > comment out warning - "No Absorption length specified"
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//              1998-08-25 by Stefano Magni
//              > Change process to use G4MaterialPropertiesTables
//              1998-09-03 by Peter Gumplinger
//              > Protect G4MaterialPropertyVector* AttenuationLengthVector
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpProcessSubType.hh"
#include "G4OpticalParameters.hh"

#include "G4OpAbsorption.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpAbsorption::G4OpAbsorption(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  Initialise();
  if(verboseLevel > 0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  SetProcessSubType(fOpAbsorption);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4OpAbsorption::~G4OpAbsorption() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpAbsorption::PreparePhysicsTable(const G4ParticleDefinition&)
{
  Initialise();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpAbsorption::Initialise()
{
  SetVerboseLevel(G4OpticalParameters::Instance()->GetAbsorptionVerboseLevel());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VParticleChange* G4OpAbsorption::PostStepDoIt(const G4Track& aTrack,
                                                const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4double thePhotonMomentum         = aParticle->GetTotalMomentum();

  aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if(verboseLevel > 1)
  {
    G4cout << "\n** OpAbsorption: Photon absorbed! **" << G4endl;
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4OpAbsorption::GetMeanFreePath(const G4Track& aTrack, G4double,
                                         G4ForceCondition*)
{
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  G4MaterialPropertiesTable* MPT =
    aTrack.GetMaterial()->GetMaterialPropertiesTable();
  G4double attLength = DBL_MAX;

  if(MPT)
  {
    G4MaterialPropertyVector* attVector = MPT->GetProperty(kABSLENGTH);
    if(attVector)
    {
      attLength =
        attVector->Value(aParticle->GetTotalMomentum(), idx_absorption);
    }
  }

  return attLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4OpAbsorption::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
  G4OpticalParameters::Instance()->SetAbsorptionVerboseLevel(verboseLevel);
}
