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

#include "G4FermiNucleiProperties.hh"

#include <G4BaryonConstructor.hh>
#include <G4NucleiProperties.hh>
#include <G4PhysicalConstants.hh>

namespace
{
// G4double WeitzsaeckerBindingEnergy(G4FermiAtomicMass atomicMass,
//   G4FermiChargeNumber chargeNumber)
// {
//   // Weitzsaecker's Mass formula
//   const auto nucleiParity = (G4int(atomicMass) - G4int(chargeNumber)) % 2;  // pairing
//   const auto chargeParity = G4int(chargeNumber) % 2;

//   const auto atomicMassF = G4double(atomicMass);
//   const auto chargeNumberF = G4double(chargeNumber);

//   G4double binding =
//     -15.67 * atomicMassF  // nuclear volume
//     + 17.23 * std::pow(atomicMassF, 2. / 3.)  // surface energy
//     + 93.15 * std::pow(atomicMassF / 2. - chargeNumberF, 2) / atomicMassF  // asymmetry
//     + 0.6984523 * std::pow(chargeNumberF, 2) / std::cbrt(atomicMassF);  // coulomb

//   if (nucleiParity == chargeParity) {
//     binding += 12.0 * (nucleiParity + chargeParity - 1) / std::sqrt(atomicMassF);  // pairing
//   }

//   return -binding * CLHEP::MeV;
// }

// G4double EstimateAtomicWeight(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
// {
//   constexpr G4double hydrogenMassExcess = 7.28897;
//   constexpr G4double neutronMassExcess = 8.07132;

//   return G4double(G4int(atomicMass) - G4int(chargeNumber)) * neutronMassExcess
//           + G4double(chargeNumber) * hydrogenMassExcess
//           - WeitzsaeckerBindingEnergy(atomicMass, chargeNumber)
//           + G4double(atomicMass) * amu_c2;
// }

// G4double EstimateNuclearMass(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
// {
//   auto mass = EstimateAtomicWeight(atomicMass, chargeNumber);

//   // atomic mass is converted to nuclear mass according formula in AME03
//   mass -= G4double(chargeNumber) * electron_mass_c2;
//   mass += (14.4381 * std::pow(G4double(chargeNumber), 2.39)
//             + 1.55468e-6 * std::pow(G4double(chargeNumber), 5.35))
//           * CLHEP::eV;

//   return mass;
// }

std::size_t GetSlot(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  const auto mass = static_cast<std::uint32_t>(atomicMass);
  const auto charge = static_cast<std::uint32_t>(chargeNumber);
  return (mass * (mass + 1)) / 2 + charge;
}
}  // namespace

G4FermiNucleiProperties::G4FermiNucleiProperties()
{
  if (G4NucleiProperties::GetNuclearMass(2, 0) <= 0.) {
    G4BaryonConstructor pCBar;
    pCBar.ConstructParticle();
  }

  for (auto a = 1; a < MAX_A; ++a) {
    for (auto z = 0; z <= a; ++z) {
      const auto atomicMass = G4FermiAtomicMass(a);
      const auto chargeNumber = G4FermiChargeNumber(z);

      const auto mass = G4NucleiProperties::GetNuclearMass(a, z);
      if (mass > 0.) {
        InsertNuclei(atomicMass, chargeNumber, mass, G4NucleiProperties::IsInStableTable(a, z));
      }
    }
  }
}

G4double G4FermiNucleiProperties::GetNuclearMassImpl(G4FermiAtomicMass atomicMass,
                                                     G4FermiChargeNumber chargeNumber) const
{
  FERMI_ASSERT_MSG(static_cast<std::uint32_t>(atomicMass)
                     >= static_cast<std::uint32_t>(chargeNumber),
                   "invalid nuclei A = " << atomicMass << ", Z = " << chargeNumber);

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (slot < nucleiMasses_.size() && nucleiMasses_[slot].isCached) {
    return nucleiMasses_[slot].mass;
  }

  return G4NucleiProperties::GetNuclearMass(G4int(atomicMass), G4int(chargeNumber));
}

G4bool G4FermiNucleiProperties::IsStableImpl(G4FermiAtomicMass atomicMass,
                                             G4FermiChargeNumber chargeNumber) const
{
  if (unlikely(atomicMass < 1_m || chargeNumber < 0_c
               || static_cast<std::uint32_t>(chargeNumber)
                    > static_cast<std::uint32_t>(atomicMass)))
  {
    return false;
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);

  return slot < nucleiMasses_.size() && nucleiMasses_[slot].isStable;
}

void G4FermiNucleiProperties::InsertNuclei(G4FermiAtomicMass atomicMass,
                                           G4FermiChargeNumber chargeNumber, G4double mass,
                                           G4bool isStable)
{
  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (slot >= nucleiMasses_.size()) {
    nucleiMasses_.resize(slot + static_cast<std::uint32_t>(atomicMass));
  }

  nucleiMasses_[slot] = G4FermiMassData{
    mass,  // mass
    isStable,  // isStable
    true,  // isCached
  };
}
