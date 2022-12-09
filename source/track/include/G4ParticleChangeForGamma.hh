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
// G4ParticleChangeForGamma
//
// Class description:
//
// Concrete class for ParticleChange for gamma processes.

// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 15 April 2005
//                                 24 August 2022
// --------------------------------------------------------------------
#ifndef G4ParticleChangeForGamma_hh
#define G4ParticleChangeForGamma_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VParticleChange.hh"

class G4DynamicParticle;

class G4ParticleChangeForGamma final : public G4VParticleChange
{
public:

  G4ParticleChangeForGamma();

  ~G4ParticleChangeForGamma() override = default;

  G4ParticleChangeForGamma(const G4ParticleChangeForGamma& right) = delete;
  G4ParticleChangeForGamma& operator=(const G4ParticleChangeForGamma& right) = delete;

  // --- the following methods are for updating G4Step -----

  G4Step* UpdateStepForAtRest(G4Step* pStep) final;
  G4Step* UpdateStepForPostStep(G4Step* Step) final;
      // A physics process gives the final state of the particle
      // based on information of G4Track

  inline void InitializeForPostStep(const G4Track&);
      // Initialize all properties by using G4Track information

  void AddSecondary(G4DynamicParticle* aParticle);
      // Add next secondary

  inline G4double GetProposedKineticEnergy() const;
  inline void SetProposedKineticEnergy(G4double proposedKinEnergy);
      // Get/Set the final kinetic energy of the current particle

  inline const G4ThreeVector& GetProposedMomentumDirection() const;
  inline void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
      // Get/Set the final momentum direction

  inline const G4ThreeVector& GetProposedPolarization() const;
  inline void ProposePolarization(const G4ThreeVector& dir);
  inline void ProposePolarization(G4double Px, G4double Py, G4double Pz);

  void DumpInfo() const override;

private:

  G4double proposedKinEnergy = 0.0;
      // The final kinetic energy of the current particle

  G4ThreeVector proposedMomentumDirection;
      // The final momentum direction of the current particle

  G4ThreeVector proposedPolarization;
      // The final polarization of the current particle
};

// ----------------------
// Inline methods
// ----------------------

inline
G4double G4ParticleChangeForGamma::GetProposedKineticEnergy() const
{
  return proposedKinEnergy;
}

inline
void G4ParticleChangeForGamma::SetProposedKineticEnergy(G4double energy)
{
  proposedKinEnergy = energy;
}

inline
const G4ThreeVector& G4ParticleChangeForGamma::
GetProposedMomentumDirection() const
{
  return proposedMomentumDirection;
}

inline
void G4ParticleChangeForGamma::
ProposeMomentumDirection(const G4ThreeVector& dir)
{
  proposedMomentumDirection = dir;
}

inline
const G4ThreeVector& G4ParticleChangeForGamma::GetProposedPolarization() const
{
  return proposedPolarization;
}

inline
void G4ParticleChangeForGamma::ProposePolarization(const G4ThreeVector& dir)
{
  proposedPolarization = dir;
}

inline
void G4ParticleChangeForGamma::ProposePolarization(G4double Px,
                                                   G4double Py,
                                                   G4double Pz)
{
  proposedPolarization.setX(Px);
  proposedPolarization.setY(Py);
  proposedPolarization.setZ(Pz);
}

inline
void G4ParticleChangeForGamma::InitializeForPostStep(const G4Track& track)
{
  InitializeSecondaries();
  InitializeLocalEnergyDeposit();
  InitializeParentWeight(track);
  InitializeStatusChange(track);
  proposedKinEnergy         = track.GetKineticEnergy();
  proposedMomentumDirection = track.GetMomentumDirection();
  proposedPolarization      = track.GetPolarization();
}

#endif
