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
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpRayleigh.hh
// Description: Discrete Process -- Rayleigh scattering of optical photons
// Version:     1.0
// Created:     1996-05-31
// Author:      Juliet Armstrong
// Updated:     2014-08-20 allow for more material types
//              2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpRayleigh_h
#define G4OpRayleigh_h 1

#include "G4VDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"

class G4OpRayleigh : public G4VDiscreteProcess
{
 public:
  explicit G4OpRayleigh(const G4String& processName = "OpRayleigh",
                        G4ProcessType type          = fOptical);
  virtual ~G4OpRayleigh();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual void BuildPhysicsTable(
    const G4ParticleDefinition& aParticleType) override;
  // Build thePhysicsTable at a right time

  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                                   G4ForceCondition*) override;
  // Returns the mean free path for Rayleigh scattering

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep) override;
  // This is the method implementing Rayleigh scattering.

  virtual G4PhysicsTable* GetPhysicsTable() const;
  // Returns the address of the physics table.

  virtual void DumpPhysicsTable() const;
  // Prints the physics table.

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;
  virtual void Initialise();

 protected:
  G4PhysicsTable* thePhysicsTable;

 private:
  G4OpRayleigh(const G4OpRayleigh& right) = delete;
  G4OpRayleigh& operator=(const G4OpRayleigh& right) = delete;

  /// Calculates the mean free paths for a material as a function of
  /// photon energy
  G4PhysicsOrderedFreeVector* CalculateRayleighMeanFreePaths(
    const G4Material* material) const;

  size_t idx_rslength = 0;
};

////////////////////
// Inline methods
////////////////////

inline G4bool G4OpRayleigh::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline void G4OpRayleigh::DumpPhysicsTable() const
{
  for(size_t i = 0; i < thePhysicsTable->entries(); ++i)
  {
    ((G4PhysicsOrderedFreeVector*) (*thePhysicsTable)[i])->DumpValues();
  }
}

inline G4PhysicsTable* G4OpRayleigh::GetPhysicsTable() const
{
  return thePhysicsTable;
}

#endif /* G4OpRayleigh_h */
