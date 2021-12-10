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
////////////////////////////////////////////////////////////////////////
// Optical Photon Absorption Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpAbsorption.hh
// Description: Discrete Process -- Bulk absorption of Optical Photons
// Version:     1.0
// Created:     1996-05-21
// Author:      Juliet Armstrong
// Updated:     2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//              1998-08-25 by Stefano Magni
//              > Change process to use G4MaterialPropertiesTables
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpAbsorption_h
#define G4OpAbsorption_h 1

#include "G4VDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"

class G4OpAbsorption : public G4VDiscreteProcess
{
 public:
  explicit G4OpAbsorption(const G4String& processName = "OpAbsorption",
                          G4ProcessType type          = fOptical);
  virtual ~G4OpAbsorption();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                                   G4ForceCondition*) override;
  // Returns the absorption length for bulk absorption of optical
  // photons in media with a specified attenuation length.

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep) override;
  // Method implementing bulk absorption of optical photons.

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;

  virtual void Initialise();

  void SetVerboseLevel(G4int);

 private:
  G4OpAbsorption(const G4OpAbsorption& right) = delete;
  G4OpAbsorption& operator=(const G4OpAbsorption& right) = delete;

  size_t idx_absorption = 0;
};

// Inline methods

inline G4bool G4OpAbsorption::IsApplicable(
  const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

#endif /* G4OpAbsorption_h */
