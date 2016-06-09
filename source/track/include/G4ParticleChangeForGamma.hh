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
// $Id: G4ParticleChangeForGamma.hh,v 1.4 2005/12/05 17:19:02 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// ------------------------------------------------------------
//	GEANT 4 class header file
//
//
// ------------------------------------------------------------
// 15 April 2005 V.Ivanchenko for gamma EM processes
//
// Modified:
// 30.05.05 : add   UpdateStepForAtRest (V.Ivanchenko)
// 04.12.05 : apply UpdateStepForPostStep in any case (mma) 
//
// ------------------------------------------------------------
//
//  Class Description
//  This class is a concrete class for ParticleChange for gamma processes
//
#ifndef G4ParticleChangeForGamma_h
#define G4ParticleChangeForGamma_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VParticleChange.hh"

class G4DynamicParticle;

class G4ParticleChangeForGamma: public G4VParticleChange
{
public:
  // default constructor
  G4ParticleChangeForGamma();

  // destructor
  virtual ~G4ParticleChangeForGamma();

  // with description
  // ----------------------------------------------------
  // --- the following methods are for updating G4Step -----

  G4Step* UpdateStepForAtRest(G4Step* pStep);
  G4Step* UpdateStepForPostStep(G4Step* Step);
  // A physics process gives the final state of the particle
  // based on information of G4Track

  void InitializeForPostStep(const G4Track&);
  //Initialize all propoerties by using G4Track information

  void AddSecondary(G4DynamicParticle* aParticle);
  // Add next secondary

  G4double GetProposedKineticEnergy() const;
  void SetProposedKineticEnergy(G4double proposedKinEnergy);
  // Get/Set the final kinetic energy of the current particle.

  const G4ThreeVector& GetProposedMomentumDirection() const;
  void ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz);
  void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
  // Get/Propose the MomentumDirection vector: it is the final momentum direction.

  const G4ThreeVector& GetProposedPolarization() const;
  void ProposePolarization(const G4ThreeVector& dir);
  void ProposePolarization(G4double Px, G4double Py, G4double Pz);

  virtual void DumpInfo() const;

  // for Debug
  virtual G4bool CheckIt(const G4Track&);

protected:
  // hide copy constructor and assignment operaor as protected
  G4ParticleChangeForGamma(const G4ParticleChangeForGamma &right);
  G4ParticleChangeForGamma & operator=(const G4ParticleChangeForGamma &right);

private:

  const G4Track* currentTrack;
  // The pointer to G4Track

  G4double proposedKinEnergy;
  //  The final kinetic energy of the current particle.

  G4ThreeVector proposedMomentumDirection;
  //  The final momentum direction of the current particle.

  G4ThreeVector proposedPolarization;
  //  The final momentum direction of the current particle.
};

// ------------------------------------------------------------

inline G4double G4ParticleChangeForGamma::GetProposedKineticEnergy() const
{
  return proposedKinEnergy;
}

inline void G4ParticleChangeForGamma::SetProposedKineticEnergy(G4double energy)
{
  proposedKinEnergy = energy;
}

inline
 const G4ThreeVector& G4ParticleChangeForGamma::GetProposedMomentumDirection() const
{
  return proposedMomentumDirection;
}

inline
 void G4ParticleChangeForGamma::ProposeMomentumDirection(const G4ThreeVector& dir)
{
  proposedMomentumDirection = dir;
}

inline
 void G4ParticleChangeForGamma::ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)
{
  proposedMomentumDirection.setX(Px);
  proposedMomentumDirection.setY(Py);
  proposedMomentumDirection.setZ(Pz);
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
 void G4ParticleChangeForGamma::ProposePolarization(G4double Px, G4double Py, G4double Pz)
{
  proposedPolarization.setX(Px);
  proposedPolarization.setY(Py);
  proposedPolarization.setZ(Pz);
}

inline void G4ParticleChangeForGamma::InitializeForPostStep(const G4Track& track)
{
  theStatusChange = track.GetTrackStatus();
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries(track);
  theParentWeight = track.GetWeight();
  proposedKinEnergy = track.GetKineticEnergy();
  proposedMomentumDirection = track.GetMomentumDirection();
  proposedPolarization = track.GetPolarization();
  currentTrack = &track;
}

//----------------------------------------------------------------
// method for updating G4Step
//

inline G4Step* G4ParticleChangeForGamma::UpdateStepForAtRest(G4Step* pStep)
{
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  pStep->SetStepLength( 0.0 );
  return pStep;
}

inline G4Step* G4ParticleChangeForGamma::UpdateStepForPostStep(G4Step* pStep)
{
  G4StepPoint* pPostStepPoint = pStep->GetPostStepPoint();
  pPostStepPoint->SetKineticEnergy( proposedKinEnergy );
  pPostStepPoint->SetMomentumDirection( proposedMomentumDirection );
  pPostStepPoint->AddPolarization( proposedPolarization );

  // update weight
  // this feature is commented out, it should be overwritten in case
  // if energy loss processes will use biasing
  // pPostStepPoint->SetWeight( theProposedWeight );
  
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
  return pStep;
}

inline void G4ParticleChangeForGamma::AddSecondary(G4DynamicParticle* aParticle)
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

