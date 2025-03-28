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

#include "G4FermiBreakUpAN.hh"

#include "G4FermiDataTypes.hh"
#include "G4FermiFragmentPool.hh"
#include "G4FermiNucleiProperties.hh"
#include "G4FermiParticle.hh"
#include "G4FermiPhaseDecay.hh"
#include "G4FermiSplitter.hh"
#include "G4FermiVFragment.hh"

#include <G4BaryonConstructor.hh>
#include <G4NucleiProperties.hh>
#include <G4PhysicalConstants.hh>
#include <G4PhysicsModelCatalog.hh>
#include <globals.hh>

#include <numeric>

#ifdef G4VERBOSE
#  define G4FERMI_VERBOSE 1
#else
#  define G4FERMI_VERBOSE 0
#endif

#define FERMI_LOG_MSG(verbosity, level, msg)                                                 \
  do {                                                                                       \
    if (G4FERMI_VERBOSE) {                                                                   \
      if ((verbosity) >= (level)) {                                                          \
        G4cout << __FILE__ << ':' << __LINE__ << " in function \"" << __FUNCTION__ << "\"\n" \
               << msg << std::endl;                                                          \
      }                                                                                      \
    }                                                                                        \
  } while (0)

constexpr G4int FERMI_WARN = 1;
constexpr G4int FERMI_DEBUG = 2;

#define FERMI_LOG_WARN(verbosity, msg) FERMI_LOG_MSG(verbosity, FERMI_WARN, msg)
#define FERMI_LOG_DEBUG(verbosity, msg) FERMI_LOG_MSG(verbosity, FERMI_DEBUG, msg)

namespace
{
constexpr const char* SPACES_OFFSET = "   ";

std::size_t SampleWeightDistribution(const std::vector<G4double>& weights)
{
  const auto totalWeight = std::accumulate(weights.begin(), weights.end(), 0.);
  FERMI_ASSERT_MSG(totalWeight > 0., "Invalid weights: all values are zero");

  const auto targetWeight = G4RandFlat::shoot() * totalWeight;
  G4double cummulativeWeight = 0;
  for (std::size_t i = 0; i < weights.size(); ++i) {
    cummulativeWeight += weights[i];

    if (cummulativeWeight >= targetWeight) {
      return i;
    }
  }

  return weights.size() - 1;
}

G4String LogProducts(const std::vector<G4FermiParticle>& particles)
{
  std::ostringstream out;

  out << "[\n";
  for (const auto& particle : particles) {
    out << SPACES_OFFSET << particle << ";\n";
  }
  out << "]";

  return std::move(out).str();
}

G4LorentzVector ChangeFrameOfReference(const G4LorentzVector& vec, const G4Vector3D& boost)
{
  auto copy = vec;
  copy.boost(boost);
  return copy;
}

G4String LogSplit(const G4FermiFragmentVector& split)
{
  std::ostringstream out;

  out << "[\n";
  for (const auto fragmentPtr : split) {
    out << SPACES_OFFSET << *fragmentPtr << ";\n";
  }
  out << "]";

  return std::move(out).str();
}

std::size_t GetSlot(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
{
  const auto mass = static_cast<std::uint32_t>(atomicMass);
  const auto charge = static_cast<std::uint32_t>(chargeNumber);
  return (mass * (mass + 1)) / 2 + charge;
}
}  // namespace

G4FermiBreakUpAN::PossibleSplits::PossibleSplits(const G4FermiAtomicMass maxAtomicMass)
{
  const auto maxMass = static_cast<std::uint32_t>(maxAtomicMass);
  splits_.resize(maxMass * (maxMass + 1) / 2);
}

const std::vector<G4FermiFragmentVector>&
G4FermiBreakUpAN::PossibleSplits::GetSplits(const G4FermiAtomicMass atomicMass,
                                            const G4FermiChargeNumber chargeNumber) const
{
  const auto slot = GetSlot(atomicMass, chargeNumber);
  return splits_.at(slot);
}

void G4FermiBreakUpAN::PossibleSplits::InsertSplits(const G4FermiAtomicMass atomicMass,
                                                    const G4FermiChargeNumber chargeNumber,
                                                    std::vector<G4FermiFragmentVector>&& splits)
{
  const auto slot = GetSlot(atomicMass, chargeNumber);

  if (slot >= splits_.size()) {
    splits_.resize(slot + static_cast<std::uint32_t>(atomicMass));
  }

  splits_[slot] = std::move(splits);
}

G4FermiBreakUpAN::G4FermiBreakUpAN(G4int verbosity)
  : splits_(G4FermiAtomicMass(MAX_A)),
    secID_(G4PhysicsModelCatalog::GetModelID("model_G4FermiBreakUpVI")),
    verbosity_(verbosity)
{}

std::vector<G4FermiParticle> G4FermiBreakUpAN::BreakItUp(const G4FermiParticle& particle) const
try {
  FERMI_LOG_DEBUG(verbosity_, "Breaking up particle: " << particle);

  if (particle.GetExcitationEnergy() < 0.) {
    FERMI_LOG_DEBUG(verbosity_, "G4FermiParticle is stable with excitation energy = "
                                  << particle.GetExcitationEnergy());
    return {particle};
  }

  const auto& splits = splits_.GetSplits(particle.GetAtomicMass(), particle.GetChargeNumber());
  FERMI_LOG_DEBUG(verbosity_,
                  "Selecting Split for " << particle << " from " << splits.size() << " splits");
  if (splits.empty()) {
    FERMI_LOG_DEBUG(verbosity_, "No splits found");
    return {particle};
  }

  // get phase space weights for every split
  // we can't cache them, because calculations is probabilistic
  weights_.resize(splits.size());
  std::transform(splits.begin(), splits.end(), weights_.begin(),
                 [atomicMass = particle.GetAtomicMass(),
                  totalEnergy = particle.GetMomentum().m()](const auto& split) {
                   return G4FermiSplitter::DecayWeight(split, atomicMass, totalEnergy);
                 });

  if (std::all_of(weights_.begin(), weights_.end(), [](auto weight) { return weight == 0.; })) {
    FERMI_LOG_DEBUG(verbosity_, "Every split has zero weight");
    return {particle};
  }

  const auto& chosenSplit = splits[SampleWeightDistribution(weights_)];
  FERMI_LOG_DEBUG(verbosity_,
                  "From " << splits.size() << " splits chosen split: " << LogSplit(chosenSplit));

  return SplitToParticles(particle, chosenSplit);
}
catch (std::exception& e) {
  G4ExceptionDescription ed;
  ed << e.what();
  G4Exception("G4FermiBreakUpAN::BreakItUp()", "Fermi002", JustWarning, ed);
  return {particle};
}

void G4FermiBreakUpAN::Initialise()
{
  if (G4NucleiProperties::GetNuclearMass(2, 0) <= 0.) {
    G4BaryonConstructor pCBar;
    pCBar.ConstructParticle();
  }
  G4FermiNucleiProperties::Instance().Initialize();

  {
    auto pool = G4FermiFragmentPool::DefaultPoolSource();
    pool.Initialize();
    G4FermiFragmentPool::Instance().Initialize(pool);
  }

  // order is important here, we use G4FermiFragmentPool to create splits!
  splits_ = PossibleSplits();
  for (auto a = 1; a < MAX_A; ++a) {
    for (auto z = 0; z <= a; ++z) {
      const auto atomicMass = G4FermiAtomicMass(a);
      const auto chargeNumber = G4FermiChargeNumber(z);

      splits_.InsertSplits(atomicMass, chargeNumber,
                           G4FermiSplitter::GenerateSplits({atomicMass, chargeNumber}));
    }
  }
}

G4bool G4FermiBreakUpAN::IsApplicable(G4int Z, G4int A, G4double /* eexc */) const
{
  return Z < MAX_Z && A < MAX_A;
}

void G4FermiBreakUpAN::BreakFragment(G4FragmentVector* results, G4Fragment* theNucleus)
{
  if (unlikely(theNucleus == nullptr || results == nullptr)) {
    G4ExceptionDescription ed;
    ed << "G4Fragment or result G4FragmentVector is not set in FermiBreakUp";
    G4Exception("G4FermiBreakUpAN::BreakFragment()", "Fermi003", FatalErrorInArgument, ed);
  }

  const auto particle =
    G4FermiParticle(G4FermiAtomicMass(theNucleus->GetA_asInt()),
                    G4FermiChargeNumber(theNucleus->GetZ_asInt()), theNucleus->GetMomentum());
  const auto fragments = BreakItUp(particle);

  const auto creationTime = theNucleus->GetCreationTime();
  for (const auto& fragment : fragments) {
    results->emplace_back(new G4Fragment(static_cast<G4int>(fragment.GetAtomicMass()),
                                         static_cast<G4int>(fragment.GetChargeNumber()),
                                         fragment.GetMomentum()));

    results->back()->SetCreationTime(creationTime);
    results->back()->SetCreatorModelID(secID_);
  }
}

std::vector<G4FermiParticle>
G4FermiBreakUpAN::SplitToParticles(const G4FermiParticle& sourceParticle,
                                   const G4FermiFragmentVector& split) const
{
  std::vector<G4double> splitMasses(split.size());
  std::transform(split.begin(), split.end(), splitMasses.begin(),
                 std::mem_fn(&G4FermiVFragment::GetTotalEnergy));

  G4FermiPhaseDecay phaseSampler;
  std::vector<G4LorentzVector> particlesMomentum;
  try {
    particlesMomentum = phaseSampler.CalculateDecay(sourceParticle.GetMomentum(), splitMasses);
  }
  catch (std::exception& e) {
    FERMI_LOG_WARN(verbosity_,
                   e.what() << " with split weight: "
                            << G4FermiSplitter::DecayWeight(split, sourceParticle.GetAtomicMass(),
                                                            sourceParticle.GetMomentum().m()));
    return {sourceParticle};
  }

  // Go back to the Lab Frame
  std::vector<G4FermiParticle> particleSplit;
  particleSplit.reserve(2 * split.size());
  const auto boostVector = sourceParticle.GetMomentum().boostVector();
  for (std::size_t fragmentIdx = 0; fragmentIdx < split.size(); ++fragmentIdx) {
    const auto fragmentMomentum =
      ChangeFrameOfReference(particlesMomentum[fragmentIdx], boostVector);
    split[fragmentIdx]->AppendDecayFragments(fragmentMomentum, particleSplit);
  }

  FERMI_LOG_DEBUG(verbosity_, "Break up products: " << LogProducts(particleSplit));
  return particleSplit;
}
