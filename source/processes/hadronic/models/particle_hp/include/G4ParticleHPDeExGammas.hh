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
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
// V.Ivanchenko 23.04.2023 Rewritten
//
#ifndef G4ParticleHPDeExGammas_h
#define G4ParticleHPDeExGammas_h 1

#include "G4ParticleHPNucLevel.hh"
#include "G4ReactionProductVector.hh"
#include "globals.hh"

#include <fstream>
#include <vector>

class G4ParticleHPDeExGammas
{
  public:
    explicit G4ParticleHPDeExGammas();
    ~G4ParticleHPDeExGammas();

    void Init(std::istream& aDataFile);

    G4ReactionProductVector* GetDecayGammas(G4int idx) const;

    inline G4int GetNumberOfLevels() const { return nLevels; }

    inline G4double GetLevelEnergy(G4int idx) const
    {
      return (idx < nLevels && idx >= 0) ? theLevels[idx]->GetLevelEnergy() : 0.0;
    }

    G4ParticleHPDeExGammas(const G4ParticleHPDeExGammas&) = delete;
    const G4ParticleHPDeExGammas& operator=(const G4ParticleHPDeExGammas&) = delete;

  private:
    G4int nLevels = 0;
    std::vector<G4ParticleHPNucLevel*> theLevels;
};

#endif
