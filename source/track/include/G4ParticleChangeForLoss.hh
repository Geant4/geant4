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
// $Id: G4ParticleChangeForLoss.hh 68795 2013-04-05 13:24:46Z gcosmo $
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
//   15.04.05 V.Ivanchenko inline update methods
//   30.01.06 V.Ivanchenko add ProposedMomentumDirection for AlongStep
//                         and ProposeWeight for PostStep
//   07.06.06 V.Ivanchenko RemoveProposedMomentumDirection from AlongStep
//   28.08.06 V.Ivanchenko Added access to current track and polarizaion
//   17.06.09 V.Ivanchenko Added SetLowEnergyLimit method 
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

  G4Step* UpdateStepForAlongStep(G4Step* Step);
  G4Step* UpdateStepForPostStep(G4Step* Step);
  // A physics process gives the final state of the particle
  // based on information of G4Track

  void InitializeForAlongStep(const G4Track&);
  void InitializeForPostStep(const G4Track&);
  //Initialize all propoerties by using G4Track information

  //  void AddSecondary(G4DynamicParticle* aParticle);
  // Add next secondary

  inline G4double GetProposedCharge() const;
  inline void SetProposedCharge(G4double theCharge);
  //   Get/Set theCharge

  inline G4double GetCharge() const;
  inline void ProposeCharge(G4double finalCharge);
  //   Get/Propose the final dynamical Charge in G4DynamicParticle

  inline G4double GetProposedKineticEnergy() const;
  inline void SetProposedKineticEnergy(G4double proposedKinEnergy);
  // Get/Set the final kinetic energy of the current particle.

  inline const G4ThreeVector& GetProposedMomentumDirection() const;
  inline void SetProposedMomentumDirection(const G4ThreeVector& dir);
  inline const G4ThreeVector& GetMomentumDirection() const;
  inline void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);
  inline void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
  // Get/Propose the MomentumDirection vector: it is the final momentum direction.

  inline const G4ThreeVector& GetProposedPolarization() const;
  inline void ProposePolarization(const G4ThreeVector& dir);
  inline void ProposePolarization(G4double Px, G4double Py, G4double Pz);

  inline const G4Track* GetCurrentTrack() const;

  inline void SetLowEnergyLimit(G4double elimit);

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

  G4double proposedKinEnergy;
  //  The final kinetic energy of the current particle.

  G4double lowEnergyLimit;
  //  The limit kinetic energy below which particle is stopped

  G4double currentCharge;
  //  The final charge of the current particle.

  G4ThreeVector proposedMomentumDirection;
  //  The final momentum direction of the current particle.

  G4ThreeVector proposedPolarization;
  //  The final polarization of the current particle.
};

// ------------------------------------------------------------

inline G4double G4ParticleChangeForLoss::GetProposedKineticEnergy() const
{
  return proposedKinEnergy;
}

inline void G4ParticleChangeForLoss::SetProposedKineticEnergy(G4double energy)
{
  proposedKinEnergy = energy;
}

inline G4double G4ParticleChangeForLoss::GetProposedCharge() const
{
  return currentCharge;
}

inline G4double G4ParticleChangeForLoss::GetCharge() const
{
  return currentCharge;
}

inline void G4ParticleChangeForLoss::SetProposedCharge(G4double theCharge)
{
  currentCharge = theCharge;
}

inline void G4ParticleChangeForLoss::ProposeCharge(G4double theCharge)
{
  currentCharge = theCharge;
}

inline
 const G4ThreeVector& G4ParticleChangeForLoss::GetProposedMomentumDirection() const
{
  return proposedMomentumDirection;
}

inline
 const G4ThreeVector& G4ParticleChangeForLoss::GetMomentumDirection() const
{
  return proposedMomentumDirection;
}

inline
 void G4ParticleChangeForLoss::ProposeMomentumDirection(const G4ThreeVector& dir)
{
  proposedMomentumDirection = dir;
}

inline
 void G4ParticleChangeForLoss::SetProposedMomentumDirection(const G4ThreeVector& dir)
{
  proposedMomentumDirection = dir;
}

inline
 void G4ParticleChangeForLoss::ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)
{
  proposedMomentumDirection.setX(Px);
  proposedMomentumDirection.setY(Py);
  proposedMomentumDirection.setZ(Pz);
}

inline const G4Track* G4ParticleChangeForLoss::GetCurrentTrack() const
{
  return currentTrack;
}

inline
 const G4ThreeVector& G4ParticleChangeForLoss::GetProposedPolarization() const
{
  return proposedPolarization;
}

inline
 void G4ParticleChangeForLoss::ProposePolarization(const G4ThreeVector& dir)
{
  proposedPolarization = dir;
}

inline
 void G4ParticleChangeForLoss::ProposePolarization(G4double Px, G4double Py, G4double Pz)
{
  proposedPolarization.setX(Px);
  proposedPolarization.setY(Py);
  proposedPolarization.setZ(Pz);
}

inline void G4ParticleChangeForLoss::InitializeForAlongStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  theNonIonizingEnergyDeposit = 0.0;
  InitializeSecondaries(track);
  theParentWeight = track.GetWeight();
  //  isParentWeightProposed = false;
  proposedKinEnergy = track.GetKineticEnergy();
  currentCharge = track.GetDynamicParticle()->GetCharge();
}

inline void G4ParticleChangeForLoss::InitializeForPostStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  theNonIonizingEnergyDeposit = 0.0;
  InitializeSecondaries(track);
  theParentWeight = track.GetWeight();
  // isParentWeightProposed = false;
  proposedKinEnergy = track.GetKineticEnergy();
  currentCharge = track.GetDynamicParticle()->GetCharge();
  proposedMomentumDirection = track.GetMomentumDirection();
  proposedPolarization = track.GetPolarization();
  currentTrack = &track;
}

/*
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
*/
inline void G4ParticleChangeForLoss::SetLowEnergyLimit(G4double elimit)
{
  lowEnergyLimit = elimit;
}

#endif

