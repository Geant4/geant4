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
// Optical Photon WaveLength Shifting (WLS) Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS2.hh
// Description: Discrete Process -- Wavelength Shifting of Optical Photons
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpWL2S_h
#  define G4OpWLS2_h 1

#  include "G4VDiscreteProcess.hh"
#  include "G4OpticalPhoton.hh"

class G4VWLSTimeGeneratorProfile;

class G4OpWLS2 : public G4VDiscreteProcess
{
 public:
  explicit G4OpWLS2(const G4String& processName = "OpWLS2",
                    G4ProcessType type          = fOptical);
  virtual ~G4OpWLS2();

  virtual G4bool IsApplicable(
    const G4ParticleDefinition& aParticleType) override;
  // Returns true -> 'is applicable' only for an optical photon.

  virtual void BuildPhysicsTable(
    const G4ParticleDefinition& aParticleType) override;
  // Build the WLS2 integral table at the right time

  virtual G4double GetMeanFreePath(const G4Track& aTrack, G4double,
                                   G4ForceCondition*) override;
  // Returns the absorption length for WLS2 absorption of optical
  // photons in media with a specified attenuation length.

  virtual G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                          const G4Step& aStep) override;
  // This is the method implementing WLS2 for optical photons.

  virtual G4PhysicsTable* GetIntegralTable() const;
  // Returns the address of the WLS2 integral table.

  virtual void DumpPhysicsTable() const;
  // Prints the WLS2 integral table.

  virtual void UseTimeProfile(const G4String name);
  // Selects the time profile generator

  virtual void PreparePhysicsTable(const G4ParticleDefinition&) override;
  virtual void Initialise();

  void SetVerboseLevel(G4int);

 protected:
  G4VWLSTimeGeneratorProfile* WLSTimeGeneratorProfile;
  G4PhysicsTable* theIntegralTable;

 private:
  G4OpWLS2(const G4OpWLS2& right) = delete;
  G4OpWLS2& operator=(const G4OpWLS2& right) = delete;

  std::size_t idx_wls2 = 0;
};

////////////////////
// Inline methods
////////////////////

inline G4bool G4OpWLS2::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return (&aParticleType == G4OpticalPhoton::OpticalPhoton());
}

inline G4PhysicsTable* G4OpWLS2::GetIntegralTable() const
{
  return theIntegralTable;
}

inline void G4OpWLS2::DumpPhysicsTable() const
{
  std::size_t PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsFreeVector* v;

  for(std::size_t i = 0; i < PhysicsTableSize; ++i)
  {
    v = (G4PhysicsFreeVector*) (*theIntegralTable)[i];
    v->DumpValues();
  }
}

#endif /* G4OpWLS2_h */
