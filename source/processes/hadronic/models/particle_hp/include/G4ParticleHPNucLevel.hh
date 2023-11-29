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
// V. Ivanchenko, 21 April 2023
//
// Data structure class for gamma levels to replace gamma data classes
// of P. Arce
//
#ifndef G4ParticleHPNucLevel_h
#define G4ParticleHPNucLevel_h 1

#include "G4ReactionProduct.hh"
#include "globals.hh"

#include <vector>

class G4ParticleHPNucLevel
{
  public:
    explicit G4ParticleHPNucLevel(G4double e);
    ~G4ParticleHPNucLevel() = default;

    void AddGamma(G4double e, G4double w, G4int idx);

    void Normalize();

    G4ReactionProduct* GetDecayGamma(G4int& idx) const;

    inline G4double GetLevelEnergy() const { return levelEnergy; }

    inline G4double GetNumberOfGammas() const { return nGammas; }

    inline G4double GetGammaEnergy(G4int idx) const
    {
      return (idx < nGammas && idx >= 0) ? gammas[idx].gammaEnergy : 0.0;
    }

    G4ParticleHPNucLevel(const G4ParticleHPNucLevel&) = delete;
    const G4ParticleHPNucLevel& operator=(const G4ParticleHPNucLevel&) = delete;

  private:
    G4int nGammas = 0;
    G4double levelEnergy;

    struct gammaData
    {
        G4double gammaEnergy;
        G4double cumProbability;
        G4int next;
    };
    std::vector<gammaData> gammas;
};

#endif
