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

#include "G4FermiBreakUp.hh"

#include "G4FermiCopyMass.hh"
#include "G4FermiDataTypes.hh"
#include "G4FermiDefaultPoolSource.hh"
#include "G4FermiFragmentPool.hh"
#include "G4FermiLogger.hh"
#include "G4FermiNucleiProperties.hh"
#include "G4FermiParticle.hh"
#include "G4FermiPhaseDecay.hh"
#include "G4FermiPossibleFragment.hh"
#include "G4FermiRandomizer.hh"
#include "G4FermiSplitter.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <algorithm>
#include <sstream>

using namespace fbu;

namespace
{
constexpr char SPACES_OFFSET[] = "   ";

G4FermiStr LogProducts(const std::vector<G4FermiParticle>& particles)
{
  std::ostringstream out;

  out << "[\n";
  for (const auto& particle : particles) {
    out << SPACES_OFFSET << particle << ";\n";
  }
  out << "]";

  return std::move(out).str();
}

G4FermiLorentzVector ChangeFrameOfReference(const G4FermiLorentzVector& vec,
                                            const G4FermiVector3& boost)
{
  auto copy = vec;
  copy.boost(boost);
  return copy;
}

std::vector<G4FermiParticle> SplitToParticles(const G4FermiParticle& sourceParticle,
                                              const G4FermiPossibleFragmentVector& split)
{
  FERMI_LOG_TRACE("Converting split to particles");

  std::vector<G4FermiFloat> splitMasses(split.size());
  std::transform(split.begin(), split.end(), splitMasses.begin(),
                 std::mem_fn(&G4FermiPossibleFragment::GetTotalEnergy));

  G4FermiPhaseDecay phaseSampler;
  std::vector<G4FermiLorentzVector> particlesMomentum;
  try {
    particlesMomentum = phaseSampler.CalculateDecay(sourceParticle.GetMomentum(), splitMasses);
  }
  catch (std::exception& e) {
    FERMI_LOG_ERROR(e.what() << " with split weight: "
                             << G4FermiSplitter::DecayWeight(split, sourceParticle.GetAtomicMass(),
                                                             sourceParticle.GetMomentum().m()));
    return {sourceParticle};
  }

  // Go back to the Lab Frame
  std::vector<G4FermiParticle> particleSplit;
  particleSplit.reserve(2 * split.size());
  const auto boostVector = sourceParticle.GetMomentum().boostVector();
  for (size_t fragmentIdx = 0; fragmentIdx < split.size(); ++fragmentIdx) {
    const auto fragmentMomentum =
      ChangeFrameOfReference(particlesMomentum[fragmentIdx], boostVector);
    split[fragmentIdx]->AppendDecayFragments(fragmentMomentum, particleSplit);
  }

  FERMI_LOG_DEBUG("Break up products: " << LogProducts(particleSplit));
  return particleSplit;
}

G4FermiStr LogSplit(const G4FermiPossibleFragmentVector& split)
{
  std::ostringstream out;

  out << "[\n";
  for (const auto fragmentPtr : split) {
    out << SPACES_OFFSET << *fragmentPtr << ";\n";
  }
  out << "]";

  return std::move(out).str();
}
}  // namespace

G4FermiBreakUp::G4FermiBreakUp(std::unique_ptr<G4FermiSplitCache>&& cache)
  : cache_(std::move(cache))
{}

std::vector<G4FermiParticle>
G4FermiBreakUp::SelectSplit(const G4FermiParticle& particle,
                            const G4FermiPossibleFragmentSplits& splits) const
{
  FERMI_LOG_TRACE("Selecting Split for " << particle << " from " << splits.size() << " splits");
  if (splits.empty()) {
    FERMI_LOG_DEBUG("No splits found");
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
    FERMI_LOG_WARN("Every split has zero weight");
    return {particle};
  }

  const auto& chosenSplit = splits[G4FermiRandomizer::SampleDistribution(weights_)];
  FERMI_LOG_DEBUG("From " << splits.size() << " splits chosen split: " << LogSplit(chosenSplit));

  return SplitToParticles(particle, chosenSplit);
}

std::vector<G4FermiParticle> G4FermiBreakUp::BreakItUp(const G4FermiParticle& particle) const
{
  FERMI_LOG_TRACE("Breaking up particle: " << particle);

  if (particle.GetExcitationEnergy() < 0.) {
    FERMI_LOG_DEBUG(
      "G4FermiParticle is stable with excitation energy = " << particle.GetExcitationEnergy());
    return {particle};
  }

  if (cache_) {
    auto splitsPtr = cache_->Get(particle.GetNucleiData());
    if (splitsPtr == nullptr) {
      FERMI_LOG_DEBUG("G4FermiParticle is stable with excitation energy = "
                      << particle.GetExcitationEnergy() / CLHEP::MeV << " MeV");
      auto splits = G4FermiSplitter::GenerateSplits(particle.GetNucleiData());
      splitsPtr = cache_->Insert(particle.GetNucleiData(), std::move(splits));
    }
    else {
      FERMI_LOG_DEBUG("Splits taken from cache");
    }

    return SelectSplit(particle, *splitsPtr);
  }
  else {
    splits_.clear();
    G4FermiSplitter::GenerateSplits(particle.GetNucleiData(), splits_);
    return SelectSplit(particle, splits_);
  }
}

void G4FermiBreakUp::Initialise()
{
  G4FermiNucleiProperties::Reset(G4FermiCopyMass());
  G4FermiFragmentPool::Reset(G4FermiDefaultPoolSource());
}

G4bool G4FermiBreakUp::IsApplicable(G4int Z, G4int A, G4double /* eexc */) const
{
  return Z < MAX_Z && A < MAX_A;
}

void G4FermiBreakUp::BreakFragment(G4FragmentVector* results, G4Fragment* theNucleus)
{
  FERMI_ASSERT_MSG(theNucleus != nullptr, "G4Fragment is not set in FermiBreakUp");
  FERMI_ASSERT_MSG(results != nullptr, "Missing result G4FragmentVector in FermiBreakUp");

  const auto particle =
    G4FermiParticle(G4FermiAtomicMass(theNucleus->GetA_asInt()),
                    G4FermiChargeNumber(theNucleus->GetZ_asInt()), theNucleus->GetMomentum());
  const auto fragments = BreakItUp(particle);

  for (const auto& fragment : fragments) {
    results->emplace_back(new G4Fragment(static_cast<G4int>(fragment.GetAtomicMass()),
                                         static_cast<G4int>(fragment.GetChargeNumber()),
                                         fragment.GetMomentum()));
  }
}
