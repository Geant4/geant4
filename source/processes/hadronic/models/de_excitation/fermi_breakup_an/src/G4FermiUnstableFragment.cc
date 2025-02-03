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
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#include "G4FermiUnstableFragment.hh"

#include "G4FermiNucleiProperties.hh"
#include "G4FermiPhaseDecay.hh"

G4FermiUnstableFragment::G4FermiUnstableFragment(G4FermiAtomicMass atomicMass,
                                                 G4FermiChargeNumber chargeNumber,
                                                 G4FermiInt polarization,
                                                 G4FermiFloat excitationEnergy,
                                                 std::vector<G4FermiNucleiData>&& decayData)
  : G4FermiVFragment(atomicMass, chargeNumber, polarization, excitationEnergy),
    decayData_(std::move(decayData))
{
  G4FermiNucleiProperties properties;
  masses_.reserve(decayData_.size());
  for (const auto& decayFragment : decayData_) {
    masses_.emplace_back(
      properties->GetNuclearMass(decayFragment.atomicMass, decayFragment.chargeNumber));
  }
}

void G4FermiUnstableFragment::AppendDecayFragments(const G4FermiLorentzVector& momentum,
                                                   std::vector<G4FermiParticle>& fragments) const
{
  G4FermiPhaseDecay phaseDecay;

  auto fragmentsMomentum = phaseDecay.CalculateDecay(momentum, masses_);

  const auto boostVector = momentum.boostVector();

  for (size_t i = 0; i < decayData_.size(); ++i) {
    fragments.emplace_back(decayData_[i].atomicMass, decayData_[i].chargeNumber,
                           fragmentsMomentum[i].boost(boostVector));
  }
}
