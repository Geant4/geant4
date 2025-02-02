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
// G4FermiBreakUp alternative de-excitation model
// by A. Novikov (January 2025)
//

#include "G4FermiFastNucleiProperties.hh"

#include "G4FermiDefaultNuclearMass.hh"
#include "G4FermiLogger.hh"

#include <CLHEP/Units/PhysicalConstants.h>

using namespace fbu;

namespace
{
G4FermiFloat WeitzsaeckerBindingEnergy(G4FermiAtomicMass atomicMass,
                                       G4FermiChargeNumber chargeNumber)
{
  // Weitzsaecker's Mass formula
  const auto nucleiParity = (G4FermiInt(atomicMass) - G4FermiInt(chargeNumber)) % 2;  // pairing
  const auto chargeParity = G4FermiInt(chargeNumber) % 2;

  const auto atomicMassF = G4FermiFloat(atomicMass);
  const auto chargeNumberF = G4FermiFloat(chargeNumber);

  G4FermiFloat binding =
    -15.67 * atomicMassF  // nuclear volume
    + 17.23 * std::pow(atomicMassF, 2. / 3.)  // surface energy
    + 93.15 * std::pow(atomicMassF / 2. - chargeNumberF, 2) / atomicMassF  // asymmetry
    + 0.6984523 * std::pow(chargeNumberF, 2) / std::cbrt(atomicMassF);  // coulomb

  if (nucleiParity == chargeParity) {
    binding += 12.0 * (nucleiParity + chargeParity - 1) / std::sqrt(atomicMassF);  // pairing
  }

  return -binding * CLHEP::MeV;
}

G4FermiFloat EstimateAtomicWeight(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  constexpr G4FermiFloat hydrogenMassExcess = 7.28897;
  constexpr G4FermiFloat neutronMassExcess = 8.07132;

  return G4FermiFloat(G4FermiInt(atomicMass) - G4FermiInt(chargeNumber)) * neutronMassExcess
         + G4FermiFloat(chargeNumber) * hydrogenMassExcess
         - WeitzsaeckerBindingEnergy(atomicMass, chargeNumber)
         + G4FermiFloat(atomicMass) * CLHEP::amu_c2;
}

G4FermiFloat EstimateNuclearMass(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  auto mass = EstimateAtomicWeight(atomicMass, chargeNumber);

  // atomic mass is converted to nuclear mass according formula in AME03
  mass -= G4FermiFloat(chargeNumber) * CLHEP::electron_mass_c2;
  mass += (14.4381 * std::pow(G4FermiFloat(chargeNumber), 2.39)
           + 1.55468e-6 * std::pow(G4FermiFloat(chargeNumber), 5.35))
          * CLHEP::eV;

  return mass;
}

inline size_t GetSlot(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  const auto mass = G4FermiUInt(atomicMass);
  const auto charge = G4FermiUInt(chargeNumber);
  return (mass * (mass + 1)) / 2 + charge;
}
}  // namespace

G4FermiFastNucleiProperties::G4FermiFastNucleiProperties()
  : G4FermiFastNucleiProperties(G4FermiDefaultNuclearMass())
{}

G4FermiFloat G4FermiFastNucleiProperties::GetNuclearMass(G4FermiAtomicMass atomicMass,
                                                         G4FermiChargeNumber chargeNumber) const
{
  FERMI_ASSERT_MSG(G4FermiUInt(atomicMass) >= G4FermiUInt(chargeNumber),
                   "invalid nuclei A = " << atomicMass << ", Z = " << chargeNumber);

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (slot < nucleiMasses_.size() && nucleiMasses_[slot].isCached) {
    return nucleiMasses_[slot].mass;
  }

  FERMI_LOG_DEBUG("Unknown particle: A = " << atomicMass << ", Z = " << chargeNumber);

  if (FERMI_UNLIKELY(slot >= nucleiMasses_.size())) {
    nucleiMasses_.resize(slot + G4FermiUInt(atomicMass));
  }

  nucleiMasses_[slot] = G4FermiMassData{
    EstimateNuclearMass(atomicMass, chargeNumber),  // mass
    false,  // isStable
    true,  // isCached
  };
  return nucleiMasses_[slot].mass;
}

bool G4FermiFastNucleiProperties::IsStable(G4FermiAtomicMass atomicMass,
                                           G4FermiChargeNumber chargeNumber) const
{
  if (FERMI_UNLIKELY(atomicMass < 1_m || chargeNumber < 0_c
                     || G4FermiUInt(chargeNumber) > G4FermiUInt(atomicMass)))
  {
    FERMI_LOG_DEBUG("Unknown particle: A = " << atomicMass << ", Z = " << chargeNumber);
    return false;
  }

  const auto slot = GetSlot(atomicMass, chargeNumber);

  return slot < nucleiMasses_.size() && nucleiMasses_[slot].isStable;
}

void G4FermiFastNucleiProperties::AddStableNuclei(G4FermiAtomicMass atomicMass,
                                                  G4FermiChargeNumber chargeNumber,
                                                  G4FermiFloat mass)
{
  FERMI_ASSERT_MSG(G4FermiUInt(atomicMass) >= G4FermiUInt(chargeNumber),
                   "invalid particle: A = " << atomicMass << ", Z = " << chargeNumber);

  const auto slot = GetSlot(atomicMass, chargeNumber);
  if (slot >= nucleiMasses_.size()) {
    nucleiMasses_.resize(slot + G4FermiUInt(atomicMass));
  }

  nucleiMasses_[slot] = G4FermiMassData{
    mass,  // mass
    true,  // isStable
    true,  // isCached
  };
}

void G4FermiFastNucleiProperties::AddStableNuclei(G4FermiNucleiData nucleiData, G4FermiFloat mass)
{
  return AddStableNuclei(nucleiData.atomicMass, nucleiData.chargeNumber, mass);
}
