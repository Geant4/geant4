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
// G4FermiBreakUpAN alternative FermiBreakUp model
// by A. Novikov (January 2025)
//
//
// Created by Artem Novikov on 30.01.2024.
//

#include "G4FermiFragmentPoolAN.hh"
#include "G4FermiDataTypes.hh"
#include "G4VFermiFragmentAN.hh"
#include "G4FermiStableFragment.hh"
#include "G4FermiUnstableFragment.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace
{
  std::size_t GetSlot(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
  {
    const auto mass = static_cast<std::uint32_t>(atomicMass);
    const auto charge = static_cast<std::uint32_t>(chargeNumber);
    return (mass * (mass + 1)) / 2 + charge;
  }
}  // namespace

G4FermiFragmentPoolAN::G4FermiFragmentPoolAN()
{
  auto pool = G4FermiFragmentPoolAN::DefaultPoolANSource();
  pool.Initialize();
  Initialize(pool);
}

std::size_t G4FermiFragmentPoolAN::Count(G4FermiAtomicMass atomicMass,
                                         G4FermiChargeNumber chargeNumber) const
{
  //  if (unlikely(static_cast<std::uint32_t>(atomicMass) < static_cast<std::uint32_t>(chargeNumber)))
  if (static_cast<std::uint32_t>(atomicMass) < static_cast<std::uint32_t>(chargeNumber)) {
    return 0;
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);
  //  if (unlikely(slot >= fragments_.size())) {
  if (slot >= fragments_.size()) {
    return 0;
  }

  return fragments_[slot].size();
}

G4FermiFragmentPoolAN::IteratorRange
G4FermiFragmentPoolAN::GetFragments(G4FermiAtomicMass atomicMass,
                                    G4FermiChargeNumber chargeNumber) const
{
  //  if (unlikely(static_cast<std::uint32_t>(atomicMass) < static_cast<std::uint32_t>(chargeNumber)))
  if (static_cast<std::uint32_t>(atomicMass) < static_cast<std::uint32_t>(chargeNumber)) {
    return {EmptyContainer_.begin(), EmptyContainer_.end()};
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (slot >= fragments_.size()) {
    return {EmptyContainer_.begin(), EmptyContainer_.end()};
  }

  return {fragments_[slot].begin(), fragments_[slot].end()};
}

void G4FermiFragmentPoolAN::AddFragment(const G4VFermiFragmentAN& fragment)
{
  const auto slot = GetSlot(fragment.GetAtomicMass(), fragment.GetChargeNumber());
  if (slot >= fragments_.size()) {
    fragments_.resize(slot + static_cast<std::uint32_t>(fragment.GetAtomicMass()));
  }
  fragments_[slot].push_back(&fragment);
}

G4FermiFragmentPoolAN::DefaultPoolANSource::DefaultPoolANSource()
{
#define FERMI_CONCAT(x, y) x##y
#define FERMI_INSTANTIATE_MACRO(x, y) FERMI_CONCAT(x, y)

#define FERMI_ADD_FRAGMENT_IMPL(NAME, VALUE) \
  static auto NAME = VALUE;                  \
  push_back(&NAME);

// automatic unique names are added
#define FERMI_ADD_FRAGMENT(VALUE) \
  FERMI_ADD_FRAGMENT_IMPL(FERMI_INSTANTIATE_MACRO(G4VFermiFragmentAN, __COUNTER__), VALUE)

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(1_m, 0_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(1_m, 1_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(2_m, 1_c, 3, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(3_m, 1_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(3_m, 2_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(4_m, 2_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(He5Fragment(5_m, 2_c, 4, 16.76 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(Li5Fragment(5_m, 3_c, 4, 16.66 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(6_m, 2_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(6_m, 3_c, 3, 0.00 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(6_m, 3_c, 1, 3.56 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(7_m, 3_c, 4, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(7_m, 3_c, 2, 0.48 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(7_m, 4_c, 4, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(7_m, 4_c, 2, 0.43 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(8_m, 3_c, 5, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(8_m, 3_c, 3, 0.98 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(Be8Fragment(8_m, 4_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(9_m, 4_c, 4, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(B9Fragment(9_m, 5_c, 4, 0.00 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 4_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 4_c, 5, 3.37 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 4_c, 8, 5.96 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 4_c, 1, 6.18 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 4_c, 5, 6.26 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 5_c, 7, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 5_c, 3, 0.72 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 5_c, 1, 1.74 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 5_c, 3, 2.15 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 5_c, 5, 3.59 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 6_c, 3, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(10_m, 6_c, 5, 3.35 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 4, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 2, 2.13 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 6, 4.44 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 4, 5.02 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 10, 6.76 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 6, 7.29 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 4, 7.98 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 5_c, 6, 8.56 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 4, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 2, 2.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 6, 4.32 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 4, 4.80 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 2, 6.34 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 8, 6.48 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 6, 6.90 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 4, 7.50 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 4, 8.10 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 6, 8.42 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(11_m, 6_c, 8, 8.66 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 5_c, 3, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 5_c, 5, 0.95 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 5_c, 5, 1.67 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 5_c, 4, 2.65 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 6_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(12_m, 6_c, 5, 4.44 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(13_m, 6_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(13_m, 6_c, 2, 3.09 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(13_m, 6_c, 4, 3.68 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(13_m, 6_c, 6, 3.85 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(13_m, 7_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 6_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 6_c, 3, 6.09 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 6_c, 8, 6.69 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 6_c, 6, 6.96 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 6_c, 5, 7.34 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 3, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 1, 2.31 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 3, 3.95 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 1, 4.92 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 5, 5.11 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 3, 5.69 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 7, 5.83 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 3, 6.20 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 7, 6.44 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(14_m, 7_c, 5, 7.03 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 8, 5.28 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 4, 6.32 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 10, 7.22 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 8, 7.57 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 2, 8.31 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 4, 8.57 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 14, 9.15 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 14, 9.79 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 7_c, 8, 10.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 8_c, 2, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 8_c, 8, 5.22 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 8_c, 4, 6.18 * CLHEP::MeV));

  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 8_c, 10, 6.83 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(15_m, 8_c, 8, 7.28 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 7_c, 5, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 7_c, 1, 0.12 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 7_c, 7, 0.30 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 7_c, 3, 0.40 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 8_c, 1, 0.00 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 8_c, 8, 6.10 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 8_c, 5, 6.92 * CLHEP::MeV));
  FERMI_ADD_FRAGMENT(G4FermiStableFragment(16_m, 8_c, 3, 7.12 * CLHEP::MeV));

#undef FERMI_ADD_FRAGMENT
#undef FERMI_ADD_FRAGMENT_IMPL
#undef FERMI_INSTANTIATE_MACRO
#undef FERMI_CONCAT
}

void G4FermiFragmentPoolAN::DefaultPoolANSource::Initialize()
{
  for (auto & fragmentPtr : *this) {
    fragmentPtr->Initialize();
  }
}
