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
// G4ParticleChangeForGamma class implementation
//
// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 15 April 2005
// --------------------------------------------------------------------

#include "G4ParticleChangeForGamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4DynamicParticle.hh"
#include "G4ExceptionSeverity.hh"

// --------------------------------------------------------------------
G4ParticleChangeForGamma::G4ParticleChangeForGamma() 
{
  // Disable flag to avoid check of each secondary at each step
  debugFlag = false;
}

// --------------------------------------------------------------------
void G4ParticleChangeForGamma::AddSecondary(G4DynamicParticle* aParticle)
{
  // create track
  G4Track* aTrack = new G4Track(aParticle, theCurrentTrack->GetGlobalTime(),
                                theCurrentTrack->GetPosition());

  // touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(theCurrentTrack->GetTouchableHandle());

  // add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForGamma::UpdateStepForAtRest(G4Step* pStep)
{
  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->SetStepLength(0.0);

  if(isParentWeightProposed)
  {
    pStep->GetPostStepPoint()->SetWeight(theParentWeight);
  }

  return pStep;
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForGamma::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  pPostStepPoint->SetMomentumDirection(proposedMomentumDirection);
  pPostStepPoint->SetPolarization(proposedPolarization);

  // update velocity for scattering process and particles with mass
  if(proposedKinEnergy > 0.0) 
  {
    pPostStepPoint->SetKineticEnergy(proposedKinEnergy);
    G4double mass = theCurrentTrack->GetDefinition()->GetPDGMass();
    G4double v = CLHEP::c_light;
    if(mass > 0.0) { 
      v *= std::sqrt(proposedKinEnergy*(proposedKinEnergy + 2*mass))/
	(proposedKinEnergy + mass);
    }
    pPostStepPoint->SetVelocity(v);
  }
  else 
  {
    pPostStepPoint->SetKineticEnergy(0.0);
    pPostStepPoint->SetVelocity(0.0);
  }

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }

  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);
  return pStep;
}

// --------------------------------------------------------------------
void G4ParticleChangeForGamma::DumpInfo() const
{
  // use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4long oldprc = G4cout.precision(8);
  G4cout << "      -----------------------------------------------" << G4endl;
  G4cout << "        G4ParticleChangeForGamma proposes: " << G4endl;
  G4cout << "        Kinetic Energy (MeV): " << std::setw(20)
         << proposedKinEnergy / MeV << G4endl;
  G4cout << "        Momentum Direction: " << std::setw(20)
         << proposedMomentumDirection << G4endl;
  G4cout << "        Polarization: " << std::setw(20) << proposedPolarization
         << G4endl;
  G4cout.precision(oldprc);
}

