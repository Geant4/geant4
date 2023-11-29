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
// G4ParticleChangeForLoss class implementation
//
// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 16 January 2004
// --------------------------------------------------------------------

#include "G4ParticleChangeForLoss.hh"
#include "G4SystemOfUnits.hh"
#include "G4ExceptionSeverity.hh"

// --------------------------------------------------------------------
G4ParticleChangeForLoss::G4ParticleChangeForLoss()
{
  // Disable flag that is enabled in G4VParticleChange if G4VERBOSE.
  debugFlag = false;
}

// --------------------------------------------------------------------
void G4ParticleChangeForLoss::DumpInfo() const
{
  // use base-class DumpInfo
  G4VParticleChange::DumpInfo();

  G4long oldprc = G4cout.precision(8);
  G4cout << "      -----------------------------------------------" << G4endl;
  G4cout << "        G4ParticleChangeForLoss proposes: " << G4endl;
  G4cout << "        Charge (eplus)   : " << std::setw(20)
         << currentCharge / eplus << G4endl;
  G4cout << "        Kinetic Energy (MeV): " << std::setw(20)
         << proposedKinEnergy / MeV << G4endl;
  G4cout << "        Momentum Direct - x : " << std::setw(20)
         << proposedMomentumDirection.x() << G4endl;
  G4cout << "        Momentum Direct - y : " << std::setw(20)
         << proposedMomentumDirection.y() << G4endl;
  G4cout << "        Momentum Direct - z : " << std::setw(20)
         << proposedMomentumDirection.z() << G4endl;
  G4cout.precision(oldprc);
}

// --------------------------------------------------------------------
G4Step* G4ParticleChangeForLoss::UpdateStepForAlongStep(G4Step* pStep)
{
  const G4StepPoint* pPreStepPoint = pStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  // accumulate change of the kinetic energy
  G4double preKinEnergy = pPreStepPoint->GetKineticEnergy();
  G4double kinEnergy =
    pPostStepPoint->GetKineticEnergy() + (proposedKinEnergy - preKinEnergy);

  pPostStepPoint->SetCharge(currentCharge);

  // calculate velocity
  if(kinEnergy > 0.0)
  {
    pPostStepPoint->SetKineticEnergy(kinEnergy);

    // assuming that mass>0, zero mass particles do not have energy loss
    pPostStepPoint->SetVelocity(CLHEP::c_light*ComputeBeta(kinEnergy));
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
G4Step* G4ParticleChangeForLoss::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();

  pPostStepPoint->SetCharge(currentCharge);
  pPostStepPoint->SetMomentumDirection(proposedMomentumDirection);
  if(proposedKinEnergy > 0.0)
  {
    pPostStepPoint->SetKineticEnergy(proposedKinEnergy);

    // assuming that mass>0, zero mass particles do not have energy loss
    pPostStepPoint->SetVelocity(CLHEP::c_light*ComputeBeta(proposedKinEnergy));
  }
  else
  {
    pPostStepPoint->SetKineticEnergy(0.0);
    pPostStepPoint->SetVelocity(0.0);
  }
  pPostStepPoint->SetPolarization(proposedPolarization);

  if(isParentWeightProposed)
  {
    pPostStepPoint->SetWeight(theParentWeight);
  }

  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);
  return pStep;
}
