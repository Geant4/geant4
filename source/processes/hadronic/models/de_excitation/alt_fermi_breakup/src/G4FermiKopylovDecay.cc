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

#include "G4FermiKopylovDecay.hh"

#include "G4FermiLogger.hh"
#include "G4FermiRandomizer.hh"

#include <CLHEP/Random/RandGamma.h>

#include <exception>
#include <functional>
#include <numeric>

using namespace fbu;

namespace
{
G4FermiLorentzVector ChangeFrameOfReference(const G4FermiLorentzVector& vec,
                                            const G4FermiVector3& boostVector)
{
  auto copy = vec;
  copy.boost(boostVector);
  return copy;
}

G4FermiFloat BetaKopylov(size_t k)
{
  constexpr G4FermiFloat beta = 1.5;

  // Notice that alpha > beta always
  const G4FermiFloat alpha = 1.5 * static_cast<G4FermiFloat>(k - 1);
  const G4FermiFloat y1 = CLHEP::RandGamma::shoot(alpha, 1);
  const G4FermiFloat y2 = CLHEP::RandGamma::shoot(beta, 1);

  return y1 / (y1 + y2);
}

G4FermiFloat TwoBodyMomentum2(G4FermiFloat totalEnergy, G4FermiFloat mass1, G4FermiFloat mass2)
{
  FERMI_ASSERT_MSG(totalEnergy > mass1 + mass2, "totalEnergy is less than fragments mass");

  return (totalEnergy + mass1 + mass2) * (totalEnergy + mass1 - mass2)
         * (totalEnergy - mass1 + mass2) * (totalEnergy - mass1 - mass2)
         / (4.0 * std::pow(totalEnergy, 2));
}

std::pair<G4FermiLorentzVector, G4FermiLorentzVector>
TwoBodyDecay(G4FermiFloat totalEnergy, G4FermiFloat mass1, G4FermiFloat mass2)
{
  const auto mag2 = TwoBodyMomentum2(totalEnergy, mass1, mass2);
  const auto momentum = G4FermiRandomizer::IsotropicVector(std::sqrt(mag2));

  return {
    G4FermiLorentzVector(momentum, std::sqrt(mag2 + std::pow(mass1, 2))),
    G4FermiLorentzVector(-momentum, std::sqrt(mag2 + std::pow(mass2, 2))),
  };
}

}  // namespace

std::vector<G4FermiLorentzVector>
G4FermiKopylovDecay::CalculateDecay(const G4FermiLorentzVector& totalMomentum,
                                    const std::vector<G4FermiFloat>& fragmentsMass) const
{
  FERMI_LOG_TRACE("Kopylov Decay called");
  FERMI_ASSERT_MSG(fragmentsMass.size() > 0, "Kopylov Decay called for empty split");

  std::vector<G4FermiLorentzVector> result(fragmentsMass.size());

  if (fragmentsMass.size() == 1) {
    FERMI_LOG_DEBUG("No decay is needed, only one fragment");
    result[0] = totalMomentum;
    return result;
  }

  // 2 bodies case is faster
  if (fragmentsMass.size() == 2) {
    FERMI_LOG_DEBUG("Decay for 2 fragments");
    std::tie(result[0], result[1]) =
      TwoBodyDecay(totalMomentum.m(), fragmentsMass[0], fragmentsMass[1]);
    return result;
  }

  // N body case
  FERMI_LOG_DEBUG("Decay for N fragments");
  auto const parentMass = totalMomentum.m();
  auto const totalFragmentsMass = std::accumulate(fragmentsMass.begin(), fragmentsMass.end(), 0.);

  auto mu = totalFragmentsMass;
  auto mass = parentMass;
  auto kineticEnergy = parentMass - totalFragmentsMass;
  FERMI_ASSERT_MSG(kineticEnergy >= 0.,
                   "Kopylov Decay started for impossible split: fragments mass is too large");

  auto momentumRestLab = G4FermiLorentzVector(0, 0, 0, parentMass);
  for (size_t i = fragmentsMass.size() - 1; i > 0; --i) {
    mu -= fragmentsMass[i];
    kineticEnergy *= i > 1 ? BetaKopylov(i) : 0.;
    const auto restMass = mu + kineticEnergy;

    FERMI_ASSERT_MSG(
      fragmentsMass[i] + restMass <= mass,
      "Kopylov Decay: something went wrong, fragments mass is greater than the whole");
    auto [momentumFragmentsCm, momentumRestCm] = TwoBodyDecay(mass, fragmentsMass[i], restMass);

    const auto boostVector = momentumRestLab.boostVector();
    momentumRestLab = ChangeFrameOfReference(momentumRestCm, boostVector);
    const auto momentumFragmentsLab = ChangeFrameOfReference(momentumFragmentsCm, boostVector);

    result[i] = momentumFragmentsLab;

    mass = restMass;
  }

  result[0] = std::move(momentumRestLab);

  return result;
}
