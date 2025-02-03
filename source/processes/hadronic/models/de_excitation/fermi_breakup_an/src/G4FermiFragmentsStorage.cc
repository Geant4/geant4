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
//
// Created by Artem Novikov on 30.01.2024.
//

#include "G4FermiFragmentsStorage.hh"

#include "G4FermiDefaultPoolSource.hh"
#include "G4FermiLogger.hh"

namespace
{
inline size_t GetSlot(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  const auto mass = G4FermiUInt(atomicMass);
  const auto charge = G4FermiUInt(chargeNumber);
  return (mass * (mass + 1)) / 2 + charge;
}
}  // namespace

G4FermiFragmentsStorage::G4FermiFragmentsStorage()
  : G4FermiFragmentsStorage(G4FermiDefaultPoolSource())
{}

size_t G4FermiFragmentsStorage::Count(G4FermiAtomicMass atomicMass,
                                      G4FermiChargeNumber chargeNumber) const
{
  if (FERMI_UNLIKELY(G4FermiUInt(atomicMass) < G4FermiUInt(chargeNumber))) {
    return 0;
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (FERMI_UNLIKELY(slot >= fragments_.size())) {
    return 0;
  }

  return fragments_[slot].size();
}

size_t G4FermiFragmentsStorage::Count(G4FermiNucleiData nuclei) const
{
  return Count(nuclei.atomicMass, nuclei.chargeNumber);
}

G4FermiFragmentsStorage::G4FermiIteratorRange
G4FermiFragmentsStorage::GetFragments(G4FermiAtomicMass atomicMass,
                                      G4FermiChargeNumber chargeNumber) const
{
  if (FERMI_UNLIKELY(G4FermiUInt(atomicMass) < G4FermiUInt(chargeNumber))) {
    return {EmptyContainer_.begin(), EmptyContainer_.end()};
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (FERMI_UNLIKELY(slot >= fragments_.size())) {
    return {EmptyContainer_.begin(), EmptyContainer_.end()};
  }

  return {fragments_[slot].begin(), fragments_[slot].end()};
}

G4FermiFragmentsStorage::G4FermiIteratorRange
G4FermiFragmentsStorage::GetFragments(G4FermiNucleiData nuclei) const
{
  return GetFragments(nuclei.atomicMass, nuclei.chargeNumber);
}

void G4FermiFragmentsStorage::AddFragment(const G4FermiVFragment& fragment)
{
  const auto slot = GetSlot(fragment.GetAtomicMass(), fragment.GetChargeNumber());
  if (slot >= fragments_.size()) {
    fragments_.resize(slot + 1);
  }
  fragments_[slot].push_back(&fragment);
}
