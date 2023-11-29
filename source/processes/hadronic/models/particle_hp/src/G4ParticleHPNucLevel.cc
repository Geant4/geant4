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
// V. Ivanchenko, 21 April 2023 Data structure class for gamma levels
//

#include "G4ParticleHPNucLevel.hh"

#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"

G4ParticleHPNucLevel::G4ParticleHPNucLevel(G4double e) : levelEnergy(e) {}

void G4ParticleHPNucLevel::AddGamma(G4double e, G4double w, G4int idx)
{
  gammaData x;
  x.gammaEnergy = e;
  x.cumProbability = w;
  x.next = idx;
  gammas.push_back(x);
  ++nGammas;
}

void G4ParticleHPNucLevel::Normalize()
{
  if (gammas.empty()) {
    return;
  }
  G4double sum = 0.0;
  for (auto& gam : gammas) {
    sum += gam.cumProbability;
  }
  if (sum <= 0.0) {
    return;
  }

  G4double norm = 1.0 / sum;
  sum = 0;
  for (auto& gam : gammas) {
    sum += norm * gam.cumProbability;
    gam.cumProbability = sum;
  }
  gammas[nGammas - 1].cumProbability = 1.0;
}

G4ReactionProduct* G4ParticleHPNucLevel::GetDecayGamma(G4int& idx) const
{
  if (gammas.empty()) {
    return nullptr;
  }
  G4double q = G4UniformRand();
  G4double e = 0.0;
  for (auto& gam : gammas) {
    if (q <= gam.cumProbability) {
      e = gam.gammaEnergy;
      idx = gam.next;
      break;
    }
  }
  if (e <= 0.0) {
    return nullptr;
  }
  G4ThreeVector p = G4RandomDirection();
  p *= e;
  auto res = new G4ReactionProduct(G4Gamma::Gamma());
  res->SetMomentum(p);
  res->SetKineticEnergy(e);
  return res;
}
