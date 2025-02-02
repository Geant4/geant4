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

#include "G4FermiSplitter.hh"

#include "G4FermiDataTypes.hh"
#include "G4FermiFragmentPool.hh"
#include "G4FermiIntegerPartition.hh"
#include "G4FermiLogger.hh"
#include "G4FermiPossibleFragment.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <optional>

using namespace fbu;

namespace
{
// Kappa = V/V_0 it is used in calculation of Coulomb energy, Kappa is dimensionless
constexpr G4FermiFloat Kappa = 1.0;

// Nuclear radius R0 (is a model parameter)
constexpr G4FermiFloat R0 = 1.3 * CLHEP::fermi;

G4FermiFloat CoulombBarrier(const G4FermiPossibleFragmentVector& split)
{
  // Coulomb Barrier (MeV) for given channel with K fragments.
  static const G4FermiFloat COEF =
    (3. / 5.) * (CLHEP::elm_coupling / R0) * std::cbrt(1. / (1. + Kappa));

  G4FermiUInt atomicMassSum = 0.;
  G4FermiUInt chargeSum = 0.;
  G4FermiFloat CoulombEnergy = 0.;
  for (const auto fragmentPtr : split) {
    auto mass = static_cast<G4FermiUInt>(fragmentPtr->GetAtomicMass());
    auto charge = static_cast<G4FermiUInt>(fragmentPtr->GetChargeNumber());
    CoulombEnergy += std::pow(charge, 2) / std::cbrt(static_cast<G4FermiFloat>(mass));
    atomicMassSum += mass;
    chargeSum += charge;
  }

  CoulombEnergy -= std::pow(static_cast<G4FermiFloat>(chargeSum), 2)
                   / std::cbrt(static_cast<G4FermiFloat>(atomicMassSum));
  return -COEF * CoulombEnergy;
}

G4FermiFloat SpinFactor(const G4FermiPossibleFragmentVector& split)
{
  G4FermiFloat factor = 1;

  for (const auto fragmentPtr : split) {
    factor *= fragmentPtr->GetPolarization();
  }

  return factor;
}

G4FermiFloat KineticEnergy(const G4FermiPossibleFragmentVector& split, G4FermiFloat totalEnergy)
{
  auto kineticEnergy = totalEnergy;
  for (const auto fragmentPtr : split) {
    kineticEnergy -= fragmentPtr->GetTotalEnergy();
  }

  // skip columb calculation for optimization purposes
  if (kineticEnergy <= 0.) {
    return kineticEnergy;
  }

  return kineticEnergy - CoulombBarrier(split);
}

G4FermiFloat MassFactor(const G4FermiPossibleFragmentVector& split)
{
  G4FermiFloat massSum = 0.;
  G4FermiFloat massProduct = 1.;
  for (const auto fragmentPtr : split) {
    const auto fragmentMass = fragmentPtr->GetMass();
    massProduct *= fragmentMass;
    massSum += fragmentMass;
  }
  auto massFactor = massProduct / massSum;
  massFactor *= std::sqrt(massFactor);
  return massFactor;
}

inline size_t Factorial(const size_t n)
{
  G4FermiUInt factorial = 1;
  for (G4FermiUInt i = 2; i <= n; ++i) {
    factorial *= i;
  }
  return factorial;
}

G4FermiFloat ConfigurationFactor(const G4FermiPossibleFragmentVector& split)
{
  // get all mass numbers and count repetitions
  std::vector<G4FermiAtomicMass> masses(split.size());
  std::transform(split.begin(), split.end(), masses.begin(),
                 std::mem_fn(&G4FermiPossibleFragment::GetAtomicMass));
  std::sort(masses.begin(), masses.end());

  // avoid overflow with floats
  // TODO: optimize with ints maybe
  G4FermiFloat factor = 1;

  size_t repeatCount = 1;  // we skip first, so start with 1
  for (size_t i = 1; i < masses.size(); ++i) {
    if (masses[i] != masses[i - 1]) {
      factor *= static_cast<G4FermiFloat>(Factorial(repeatCount));
      repeatCount = 0;
    }
    ++repeatCount;
  }
  factor *= static_cast<G4FermiFloat>(Factorial(repeatCount));

  return factor;
}

G4FermiFloat ConstFactor(G4FermiAtomicMass atomicMass, size_t fragmentsCount)
{
  static const G4FermiFloat COEF =
    std::pow(R0 / CLHEP::hbarc, 3) * Kappa * std::sqrt(2.0 / CLHEP::pi) / 3.0;

  return std::pow(COEF * static_cast<G4FermiFloat>(atomicMass), fragmentsCount - 1);
}

G4FermiFloat GammaFactor(size_t fragmentsCount)
{
  G4FermiFloat gamma = 1.0;
  G4FermiFloat arg = 3.0 * static_cast<G4FermiFloat>(fragmentsCount - 1) / 2.0 - 1.0;
  while (arg > 1.1) {
    gamma *= arg;
    arg -= 1;
  }

  if (fragmentsCount % 2 == 0) {
    gamma *= std::sqrt(CLHEP::pi);
  }

  return gamma;
}
}  // namespace

G4FermiFloat G4FermiSplitter::DecayWeight(const G4FermiPossibleFragmentVector& split,
                                          G4FermiAtomicMass atomicMass, G4FermiFloat totalEnergy)
{
  const auto kineticEnergy = KineticEnergy(split, totalEnergy);  // in MeV

  // Check that there is enough energy to produce K fragments
  if (kineticEnergy <= 0.) {
    return 0.;
  }

  const auto power = 3.0 * static_cast<G4FermiFloat>(split.size() - 1) / 2.0 - 1.;
  const auto kineticFactor = std::pow(kineticEnergy, power);

  // Spin factor S_n
  const auto spinFactor = SpinFactor(split);

  // Calculate MassFactor
  const auto massFactor = MassFactor(split);

  // This is the constant (doesn't depend on energy) part
  const auto coef = ConstFactor(atomicMass, split.size());

  // Calculation of 1/gamma(3(k-1)/2)
  const auto gamma = GammaFactor(split.size());

  // Permutation Factor G_n
  const auto permutationFactor = ConfigurationFactor(split);

  return coef * kineticFactor * massFactor * spinFactor / (permutationFactor * gamma);
}

namespace
{
constexpr size_t ExpectedSplitSize = 100;

void ThrowOnInvalidInputs(G4FermiNucleiData nucleiData)
{
  FERMI_ASSERT_MSG(nucleiData.atomicMass > 0_m && nucleiData.chargeNumber >= 0_c,
                   "Non valid arguments A = " << nucleiData.atomicMass
                                              << " Z = " << nucleiData.chargeNumber);

  FERMI_ASSERT_MSG(G4FermiUInt(nucleiData.chargeNumber) <= G4FermiUInt(nucleiData.atomicMass),
                   "Non physical arguments = " << nucleiData.atomicMass
                                               << " Z = " << nucleiData.chargeNumber);
}

G4FermiPossibleFragmentSplits PossibleSplits(const G4FermiPartition& massPartition,
                                             const G4FermiPartition& chargePartition)
{
  auto& fragmentPool = G4FermiFragmentPool::Instance();
  const auto fragmentCount = massPartition.size();

  // count all possible splits due to multiplicity of fragments
  size_t splitsCount = 1;
  for (size_t fragmentIdx = 0; fragmentIdx < fragmentCount; ++fragmentIdx) {
    splitsCount *= fragmentPool.Count(G4FermiAtomicMass(massPartition[fragmentIdx]),
                                      G4FermiChargeNumber(chargePartition[fragmentIdx]));
    if (splitsCount == 0) {
      return {};
    }
  }

  // allocate in advance
  G4FermiPossibleFragmentSplits splits(splitsCount, G4FermiPossibleFragmentVector(fragmentCount));

  // incrementally build splits
  // !! chosen order matters, because later there we need to remove duplicates
  size_t groupSize = splitsCount;
  for (size_t fragmentIdx = 0; fragmentIdx < fragmentCount; ++fragmentIdx) {
    const auto fragmentRange =
      fragmentPool.GetFragments(G4FermiAtomicMass(massPartition[fragmentIdx]),
                                G4FermiChargeNumber(chargePartition[fragmentIdx]));
    // no remainder here!
    const size_t multiplicity = std::distance(fragmentRange.begin(), fragmentRange.end());
    groupSize /= multiplicity;

    for (size_t offset = 0; offset < splitsCount;) {
      for (const auto fragmentPtr : fragmentRange) {
        for (size_t pos = 0; pos < groupSize; ++pos) {
          splits[offset + pos][fragmentIdx] = fragmentPtr;
        }
        offset += groupSize;
      }
    }
  }

  // remove duplicate splits
  for (auto& split : splits) {
    std::sort(split.begin(), split.end(), std::greater<>());
    // greater, because they already partially sorted as greater due to integer partition
  }
  const auto uniqueEndIt = std::unique(splits.begin(), splits.end());
  splits.resize(uniqueEndIt - splits.begin());

  return splits;
}
}  // namespace

G4FermiPossibleFragmentSplits G4FermiSplitter::GenerateSplits(G4FermiNucleiData nucleiData)
{
  G4FermiPossibleFragmentSplits splits;
  GenerateSplits(nucleiData, splits);
  return splits;
}

void G4FermiSplitter::GenerateSplits(G4FermiNucleiData nucleiData,
                                     G4FermiPossibleFragmentSplits& splits)
{
  ThrowOnInvalidInputs(nucleiData);

  splits.reserve(ExpectedSplitSize);

  // let's split nucleus into 2, ..., A fragments
  auto maxFragmentsCount = G4FermiUInt(nucleiData.atomicMass);

  for (G4FermiUInt fragmentCount = 2; fragmentCount <= maxFragmentsCount; ++fragmentCount) {
    // Form all possible partition by combination of A partitions and Z partitions (Z partitions
    // include null parts)
    for (auto& massPartition : G4FermiIntegerPartition(nucleiData.atomicMass, fragmentCount, 1)) {
      for (auto& chargePartition :
           G4FermiIntegerPartition(nucleiData.chargeNumber, fragmentCount, 0)) {
        // Some splits are invalid, some nuclei doesn't exist
        if (auto partitionSplits = PossibleSplits(massPartition, chargePartition);
            !partitionSplits.empty()) {
          splits.insert(splits.end(), std::make_move_iterator(partitionSplits.begin()),
                        std::make_move_iterator(partitionSplits.end()));
        }
      }
    }
  }
}
