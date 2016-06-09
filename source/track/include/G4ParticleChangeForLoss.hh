//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ParticleChangeForLoss.hh,v 1.5 2004/01/20 15:29:41 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-01 $
//
//
// ------------------------------------------------------------
//	GEANT 4 class header file
//
//
// ------------------------------------------------------------
//   Implemented for the new scheme                 23 Mar. 1998  H.Kurahige
//
//   Modified:
//   16.01.04 V.Ivanchenko update for model variant of energy loss
//
// ------------------------------------------------------------
//
//  Class Description
//  This class is a concrete class for ParticleChange for EnergyLoss
//
#ifndef G4ParticleChangeForLoss_h
#define G4ParticleChangeForLoss_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VParticleChange.hh"

class G4DynamicParticle;

class G4ParticleChangeForLoss: public G4VParticleChange
{
public:
  // default constructor
  G4ParticleChangeForLoss();

  // destructor
  virtual ~G4ParticleChangeForLoss();

  // with description
  // ----------------------------------------------------
  // --- the following methods are for updating G4Step -----

  virtual G4Step* UpdateStepForAlongStep(G4Step* Step);
  virtual G4Step* UpdateStepForPostStep(G4Step* Step);
  // A physics process gives the final state of the particle
  // based on information of G4Track

  void InitializeForAlongStep(const G4Track&);
  void InitializeForPostStep(const G4Track&);
  //Initialize all propoerties by using G4Track information

  void AddSecondary(G4DynamicParticle* aParticle);
  // Add next secondary

  G4double GetProposedCharge() const;
  void SetProposedCharge(G4double theCharge);
  //   Get/Set theCharge

  G4double GetProposedKineticEnergy() const;
  void SetProposedKineticEnergy(G4double kinEnergy);
  // Get/Set the final kinetic energy of the current particle.

  const G4ThreeVector& GetProposedMomentumDirection() const;
  void SetProposedMomentumDirection(const G4ThreeVector& dir);

  virtual void DumpInfo() const;

  // for Debug
  virtual G4bool CheckIt(const G4Track&);

protected:
  // hide copy constructor and assignment operaor as protected
  G4ParticleChangeForLoss(const G4ParticleChangeForLoss &right);
  G4ParticleChangeForLoss & operator=(const G4ParticleChangeForLoss &right);

private:

  const G4Track* currentTrack;
  // The pointer to G4Track
  G4double kinEnergy;
  //  The final kinetic energy of the current particle.
  G4double currentCharge;
  //  The final charge of the current particle.
  G4ThreeVector proposedMomentumDirection;
  //  The final momentum direction of the current particle.
};

inline G4double G4ParticleChangeForLoss::GetProposedKineticEnergy() const
{
  return kinEnergy;
}

inline void G4ParticleChangeForLoss::SetProposedKineticEnergy(G4double energy)
{
  kinEnergy = energy;
}

inline G4double G4ParticleChangeForLoss::GetProposedCharge() const
{
  return currentCharge;
}

inline void G4ParticleChangeForLoss::SetProposedCharge(G4double theCharge)
{
  currentCharge = theCharge;
}

inline
 const G4ThreeVector& G4ParticleChangeForLoss::GetProposedMomentumDirection() const
{
  return proposedMomentumDirection;
}

inline
 void G4ParticleChangeForLoss::SetProposedMomentumDirection(const G4ThreeVector& dir)
{
  proposedMomentumDirection = dir;
}

inline void G4ParticleChangeForLoss::InitializeForAlongStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries(track);
  kinEnergy = track.GetKineticEnergy();
  currentCharge = track.GetDynamicParticle()->GetCharge();
}

inline void G4ParticleChangeForLoss::InitializeForPostStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries(track);
  kinEnergy = track.GetKineticEnergy();
  currentCharge = track.GetDynamicParticle()->GetCharge();
  proposedMomentumDirection = track.GetMomentumDirection();
  currentTrack = &track;
}

inline void G4ParticleChangeForLoss::AddSecondary(G4DynamicParticle* aParticle)
{
  //  create track
  G4Track* aTrack = new G4Track(aParticle, currentTrack->GetGlobalTime(),
                                           currentTrack->GetPosition());

  //   Touchable handle is copied to keep the pointer
  aTrack->SetTouchableHandle(currentTrack->GetTouchableHandle());

  //  add a secondary
  G4VParticleChange::AddSecondary(aTrack);
}

#endif

