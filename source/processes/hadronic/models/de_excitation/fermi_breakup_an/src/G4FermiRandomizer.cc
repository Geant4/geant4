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

#include "G4FermiRandomizer.hh"

#include "G4FermiLogger.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <algorithm>
#include <numeric>
#include <random>

G4FermiFloat G4FermiRandomizer::SampleUniform()
{
  static std::uniform_real_distribution<G4FermiFloat> distribution(0., 1.);
  return distribution(Engine_);
}

G4FermiFloat G4FermiRandomizer::SampleNormal(G4FermiFloat mean, G4FermiFloat deviation)
{
  static std::normal_distribution<G4FermiFloat> distribution(0., 1.);
  return mean + distribution(Engine_) * deviation;
}

G4FermiParticleMomentum G4FermiRandomizer::IsotropicVector(G4FermiFloat magnitude)
{
  const auto cos = 1.0 - 2.0 * SampleUniform();
  const auto sin = std::sqrt(1.0 - std::pow(cos, 2));
  const auto phi = CLHEP::twopi * SampleUniform();

  return G4FermiParticleMomentum(magnitude * std::cos(phi) * sin, magnitude * std::sin(phi) * sin,
                                 magnitude * cos);
}

std::vector<G4FermiFloat> G4FermiRandomizer::ProbabilityDistribution(size_t pointCount)
{
  // Sample uniform random numbers in increasing order
  std::vector<G4FermiFloat> probabilityDistribution(pointCount);

  probabilityDistribution.front() = 0.;
  std::generate_n(probabilityDistribution.begin(), pointCount - 2,
                  G4FermiRandomizer::SampleUniform);
  probabilityDistribution.back() = 1.;

  std::sort(probabilityDistribution.begin(), probabilityDistribution.end());

  return probabilityDistribution;
}

size_t G4FermiRandomizer::SampleDistribution(const std::vector<G4FermiFloat>& weights)
{
  const auto totalWeight = std::accumulate(weights.begin(), weights.end(), 0.);
  FERMI_ASSERT_MSG(totalWeight > 0., "Invalid weights: all values are zero");

  const auto targetWeight = SampleUniform() * totalWeight;
  G4FermiFloat cummulativeWeight = 0;
  for (size_t i = 0; i < weights.size(); ++i) {
    cummulativeWeight += weights[i];

    if (cummulativeWeight >= targetWeight) {
      return i;
    }
  }

  return weights.size() - 1;
}

void G4FermiRandomizer::SetSeed(G4FermiRandomEngine::result_type seed)
{
  Engine_.seed(seed);
}
