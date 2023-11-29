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
// G4ParticleChangeForMSC
//
// Class description:
//
// Concrete "Particle Change" class for Multiple Scattering process.

// Author: Hisaya Kurashige, 23 March 1998  
// Revision: Vladimir Ivantchenko, 16 January 2004
//                                 23 August 2022
// --------------------------------------------------------------------
#ifndef G4ParticleChangeForMSC_hh
#define G4ParticleChangeForMSC_hh 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include "G4VParticleChange.hh"

class G4ParticleChangeForMSC final : public G4VParticleChange
{
public:

  G4ParticleChangeForMSC();

  ~G4ParticleChangeForMSC() override = default;

  G4ParticleChangeForMSC(const G4ParticleChangeForMSC& right) = delete;
  G4ParticleChangeForMSC& operator=
  (const G4ParticleChangeForMSC& right) = delete;

  // A multiple scattering process gives the final state of the particle
  G4Step* UpdateStepForAlongStep(G4Step* step) final;

  // Initialise only properties changed by multiple scattering
  inline void InitialiseMSC(const G4Track&, const G4Step& step);

  inline void ProposeMomentumDirection(const G4ThreeVector& Pfinal);
  inline const G4ThreeVector* GetProposedMomentumDirection() const;

  inline void ProposePosition(const G4ThreeVector& finalPosition);
  inline const G4ThreeVector& GetProposedPosition() const;

private:

  G4ThreeVector theMomentumDirection;
  // It is the vector containing the final momentum direction
  // after multiple scattering is applied along step

  G4ThreeVector thePosition;
  // The changed of post step position after multiple scattering
};

inline void G4ParticleChangeForMSC::ProposeMomentumDirection(
  const G4ThreeVector& P)
{
  theMomentumDirection = P;
}

inline const G4ThreeVector*
G4ParticleChangeForMSC::GetProposedMomentumDirection() const
{
  return &theMomentumDirection;
}

inline void G4ParticleChangeForMSC::ProposePosition(const G4ThreeVector& P)
{
  thePosition = P;
}

inline const G4ThreeVector& G4ParticleChangeForMSC::GetProposedPosition() const
{
  return thePosition;
}

inline void
G4ParticleChangeForMSC::InitialiseMSC(const G4Track& track, const G4Step& step)
{
  theStatusChange = track.GetTrackStatus();
  auto poststep = step.GetPostStepPoint();
  thePosition = poststep->GetPosition();
  theMomentumDirection = poststep->GetMomentumDirection();  
  theCurrentTrack = &track;
}

#endif
