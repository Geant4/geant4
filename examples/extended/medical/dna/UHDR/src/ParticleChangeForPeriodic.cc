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
/*
 * Based on 'G4pbc'.
 * Copyright (c) 2020 Amentum Pty Ltd
 * team@amentum.space
 * The original open-source version of this code
 * may be found at https://github.com/amentumspace/g4pbc
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 * associated documentation files (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge, publish, distribute,
 * sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
 * is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies
 * or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 * NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 */
#include "ParticleChangeForPeriodic.hh"

#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ParticleChangeForPeriodic::ParticleChangeForPeriodic() : G4VParticleChange() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Step* ParticleChangeForPeriodic::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  pPostStepPoint->SetMomentumDirection(fProposedMomentumDirection);
  pPostStepPoint->SetPolarization(fProposedPolarization);
  pPostStepPoint->SetPosition(fProposedPosition);

  if (isParentWeightProposed) {
    pPostStepPoint->SetWeight(theParentWeight);
  }

  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);

  return pStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticleChangeForPeriodic::AddSecondary(G4DynamicParticle* aParticle)
{
  auto* aTrack = new G4Track(aParticle, fTrack->GetGlobalTime(), fTrack->GetPosition());

  aTrack->SetTouchableHandle(fTrack->GetTouchableHandle());

  G4VParticleChange::AddSecondary(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticleChangeForPeriodic::DumpInfo() const
{
  G4VParticleChange::DumpInfo();
  G4int oldprc = G4cout.precision(3);

  G4cout << "        Momentum Direction: " << std::setw(20) << fProposedMomentumDirection << G4endl;
  G4cout << "        Polarization: " << std::setw(20) << fProposedPolarization << G4endl;
  G4cout << "        Position: " << std::setw(20) << fProposedPosition << G4endl;
  G4cout.precision(oldprc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
