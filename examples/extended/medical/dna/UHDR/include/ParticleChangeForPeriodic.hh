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
#ifndef ParticleChangeForPeriodic_hh
#define ParticleChangeForPeriodic_hh 1

#include "G4VParticleChange.hh"
#include "globals.hh"

class G4DynamicParticle;

class ParticleChangeForPeriodic : public G4VParticleChange
{
  public:
    ParticleChangeForPeriodic();

    ~ParticleChangeForPeriodic() override = default;

    G4Step* UpdateStepForPostStep(G4Step* Step) override;

    void InitializeForPostStep(const G4Track&);

    void AddSecondary(G4DynamicParticle* aParticle);

    const G4ThreeVector& GetProposedMomentumDirection() const;

    void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);

    void ProposeMomentumDirection(const G4ThreeVector& Pfinal);

    const G4ThreeVector& GetProposedPolarization() const;

    void ProposePolarization(const G4ThreeVector& dir);

    void ProposePolarization(G4double Px, G4double Py, G4double Pz);

    const G4ThreeVector& GetProposedPosition() const;

    void ProposePosition(const G4ThreeVector& pos);

    void ProposePosition(G4double x, G4double y, G4double z);

    const G4Track* GetCurrentTrack() const;

    void DumpInfo() const override;

    ParticleChangeForPeriodic(const ParticleChangeForPeriodic& right) = delete;

    ParticleChangeForPeriodic& operator=(const ParticleChangeForPeriodic& right) = delete;

  private:
    const G4Track* fTrack{};
    G4ThreeVector fProposedMomentumDirection;
    G4ThreeVector fProposedPolarization;
    G4ThreeVector fProposedPosition;
};

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedMomentumDirection() const
{
  return fProposedMomentumDirection;
}

inline void ParticleChangeForPeriodic::ProposeMomentumDirection(const G4ThreeVector& dir)
{
  fProposedMomentumDirection = dir;
}

inline void ParticleChangeForPeriodic::ProposeMomentumDirection(G4double Px, G4double Py,
                                                                G4double Pz)
{
  fProposedMomentumDirection.setX(Px);
  fProposedMomentumDirection.setY(Py);
  fProposedMomentumDirection.setZ(Pz);
}

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedPolarization() const
{
  return fProposedPolarization;
}

inline void ParticleChangeForPeriodic::ProposePolarization(const G4ThreeVector& dir)
{
  fProposedPolarization = dir;
}

inline void ParticleChangeForPeriodic::ProposePolarization(G4double Px, G4double Py, G4double Pz)
{
  fProposedPolarization.setX(Px);
  fProposedPolarization.setY(Py);
  fProposedPolarization.setZ(Pz);
}

inline const G4ThreeVector& ParticleChangeForPeriodic::GetProposedPosition() const
{
  return fProposedPosition;
}

inline void ParticleChangeForPeriodic::ProposePosition(const G4ThreeVector& dir)
{
  fProposedPosition = dir;
}

inline void ParticleChangeForPeriodic::ProposePosition(G4double Px, G4double Py, G4double Pz)
{
  fProposedPosition.setX(Px);
  fProposedPosition.setY(Py);
  fProposedPosition.setZ(Pz);
}

inline void ParticleChangeForPeriodic::InitializeForPostStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  theNonIonizingEnergyDeposit = 0.0;
  InitializeSecondaries();
  theParentWeight = track.GetWeight();
  isParentWeightProposed = false;
  fProposedMomentumDirection = track.GetMomentumDirection();
  fProposedPolarization = track.GetPolarization();
  fProposedPosition = track.GetPosition();
  fTrack = &track;
}

inline const G4Track* ParticleChangeForPeriodic::GetCurrentTrack() const
{
  return fTrack;
}
#endif