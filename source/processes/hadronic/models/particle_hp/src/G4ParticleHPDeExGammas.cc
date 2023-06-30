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
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
// 080801 Prohibit level transition to oneself by T. Koi
//
// P. Arce, June-2014 Conversion neutron_hp to particle_hp
//
#include "G4ParticleHPDeExGammas.hh"

#include "G4RandomDirection.hh"
#include "G4ReactionProduct.hh"
#include "G4SystemOfUnits.hh"

G4ParticleHPDeExGammas::G4ParticleHPDeExGammas() = default;

G4ParticleHPDeExGammas::~G4ParticleHPDeExGammas()
{
  for (auto& ptr : theLevels) {
    delete ptr;
  }
}

void G4ParticleHPDeExGammas::Init(std::istream& aDataFile)
{
  // G4cout << "### G4ParticleHPDeExGammas::Init new file " <<  G4endl;
  // ground state
  auto level = new G4ParticleHPNucLevel(0.0);

  G4double elevel0 = 0.0;
  G4double elevel = 0.0;
  G4double egamma = 0.0;
  G4double prob = 0.0;
  constexpr G4double eps = 1 * CLHEP::eV;
  for (;;) {
    if (aDataFile >> elevel) {
      // next line
      aDataFile >> egamma >> prob;
      egamma *= CLHEP::keV;
      elevel *= CLHEP::keV;
      prob = std::max(prob, 1.e-6);
      // G4cout << "    El0=" << elevel0 << " El=" << elevel
      //  << " Eg=" << egamma << " w=" << prob << G4endl;

      // save previous level and start a new level
      if (std::abs(elevel - elevel0) > eps) {
        level->Normalize();
        theLevels.push_back(level);
        ++nLevels;
        level = new G4ParticleHPNucLevel(elevel);
        elevel0 = elevel;
        // G4cout << "  New level " << nLevels << "  E=" << elevel << G4endl;
      }

      // find the next level
      G4double e = elevel - egamma;
      G4int next = -1;
      G4double del = DBL_MAX;
      for (G4int i = 0; i < nLevels; ++i) {
        G4double de = std::abs(theLevels[i]->GetLevelEnergy() - e);
        if (de < del) {
          next = i;
          del = de;
        }
      }
      // save level data
      if (next >= 0) {
        level->AddGamma(egamma, prob, next);
        // G4cout << "      NLevel=" << nLevels << " Elevel=" << elevel
        //      << " Egamma=" << egamma << " W=" << prob << " next=" << next << G4endl;
      }
    }
    else {
      // end of file - save recent level
      level->Normalize();
      theLevels.push_back(level);
      ++nLevels;
      // G4cout << "### End of file Nlevels=" << nLevels << G4endl;
      break;
    }
  }
}

G4ReactionProductVector* G4ParticleHPDeExGammas::GetDecayGammas(G4int i) const
{
  G4int idx = i;
  if (idx >= nLevels || idx <= 0) return nullptr;
  auto result = new G4ReactionProductVector();

  for (;;) {
    if (idx <= 0) {
      break;
    }
    auto ptr = theLevels[idx]->GetDecayGamma(idx);
    if (nullptr != ptr) {
      result->push_back(ptr);
    }
  }
  return result;
}
