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
// G4MuonDecayChannel class implementation
//
// Author: H.Kurashige, 30 May 1997
// Contributions:
// - 2005 - M.Melissas, J.Brunner CPPM/IN2P3
//   Added V-A fluxes for neutrinos using a new algorithm, 2005
// - S.Okada KEK/CRC, 14 January 2026
//   Fixed energy conservation in G4MuonDecayChannel (Bugzilla #2678):
//   Corrected EMax definition to parentmass/2 and applied standard relativistic
//   formula (p=sqrt(E^2-m^2)) to compute momenta for all decay products.
// --------------------------------------------------------------------

#include "G4MuonDecayChannel.hh"

#include "G4DecayProducts.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "G4VDecayChannel.hh"
#include "Randomize.hh"

#include <algorithm>

namespace
{
/**
 * \brief Compute \f$\tan(C/2)\f$ for the interior angle \f$C\f$ opposite side \f$c\f$.
 *
 * For a triangle with side lengths \f$(a,b,c)\f$, this helper evaluates a numerically stable
 * form of the tangent half-angle relation for \f$C\f$:
 *
 * \f[
 * \tan\left(\frac{C}{2}\right)
 * = \sqrt{\frac{((a-b)+c)\,\mu}{(a+(b+c))\,((a-c)+b)}}
 * \f]
 *
 * where
 * \f[
 * \mu =
 * \begin{cases}
 * c-(a-b), & b \ge c \\
 * b-(a-c), & b < c
 * \end{cases}
 * \f]
 *
 * and \f$a\f$ and \f$b\f$ are internally reordered so \f$a \ge b\f$ before evaluation.
 * This arrangement mitigates cancellation for near-degenerate triangles.
 *
 * Reference:
 * - W. Kahan, "Miscalculating Area and Angles of a Needle-like Triangle"
 *   https://people.eecs.berkeley.edu/~wkahan/Triangle.pdf
 *
 * \pre \f$(a,b,c)\f$ must satisfy triangle closure (all side lengths positive and
 *      \f$|a-b| \le c \le a+b\f$).
 * \param a Triangle side length.
 * \param b Triangle side length.
 * \param c Triangle side length opposite angle \f$C\f$.
 * \return \f$\tan(C/2)\f$.
 */
inline G4double kahan_interior_angle(G4double a, G4double b, G4double c)
{
  if (a < b)
  {
    std::swap(a, b);
  }
  // Parentheses _must_ be retained as these are part of the numerical
  // stability.
  const G4double mu = (b >= c) ? c - (a - b) : b - (a - c);
  const G4double numerator = ((a - b) + c) * mu;
  const G4double denominator = (a + (b + c)) * ((a - c) + b);
  return std::sqrt(std::max(0.0, numerator / denominator));
}
}  // namespace

G4MuonDecayChannel::G4MuonDecayChannel(const G4String& theParentName, G4double theBR)
  : G4VDecayChannel("Muon Decay", 1)
{
  // set names for daughter particles
  if (theParentName == "mu+")
  {
    SetBR(theBR);
    SetParent("mu+");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e+");
    SetDaughter(1, "nu_e");
    SetDaughter(2, "anti_nu_mu");
  }
  else if (theParentName == "mu-")
  {
    SetBR(theBR);
    SetParent("mu-");
    SetNumberOfDaughters(3);
    SetDaughter(0, "e-");
    SetDaughter(1, "anti_nu_e");
    SetDaughter(2, "nu_mu");
  }
  else
  {
#ifdef G4VERBOSE
    if (GetVerboseLevel() > 0)
    {
      G4cout << "G4MuonDecayChannel:: constructor :";
      G4cout << " parent particle is not muon but ";
      G4cout << theParentName << G4endl;
    }
#endif
  }
}

G4MuonDecayChannel& G4MuonDecayChannel::operator=(const G4MuonDecayChannel& right)
{
  if (this != &right)
  {
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if (numberOfDaughters > 0)
    {
      if (daughters_name != nullptr) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      // copy daughters name
      for (G4int index = 0; index < numberOfDaughters; ++index)
      {
        daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
  }
  return *this;
}

/**
 * \brief Samples the three-body decay \f$ \mu^- \rightarrow e^- + \bar{\nu_e} + \nu_{\mu} \f$ and
 * its antiparticle equivalent.
 * \ingroup particles_management
 *
 * This implementation assumes the Standard Model V-A matrix element and phase space, neglecting
 * muon polarization and neutrino masses, but retaining the electron mass. Here,
 * the spin-summed matrix element in terms of the particle four-momenta and corresponding
 * differential decay rate are:
 *
 * \f[
 * |\mathcal{M}|^2 = 64 G_F^2 (p_\mu \cdot p_{\bar{\nu}_e})(p_e \cdot p_{\nu_\mu}),
 * \qquad
 * d\Gamma = \frac{1}{2m_\mu} |\mathcal{M}|^2 d\Phi_3,
 * \f]
 *
 * where \f$d\Phi_3\f$ is the Lorentz-invariant three-body phase-space element.
 * Taking reduced energies as:
 *
 *\f[
 * x = \frac{2E_e}{m_\mu},
 * \qquad
 * y = \frac{2E_{\bar{\nu_e}}}{m_\mu},
 * \f]
 *
 * it can be shown that integrating the double-differential decay rate over the
 * electron energy yields the marginal electron-antineutrino spectrum
 *
 * \f[
 * \frac{d\Gamma}{dy}
 * \propto
 * \frac{y^2(a-y)^2}{1-y},
 * \f]
 *
 * where
 *
 * \f[
 * r = \frac{m_e}{m_\mu},
 * \qquad
 * a = 1-r^2.
 * \f]
 *
 * The electron-antineutrino energy is sampled using acceptance-rejection from
 * this distribution. The rejection envelope is constant over \f$y\f$ with a value
 * of \f$27/4\f$, the maximum of \f$d\Gamma/dy\f$ when the electron mass is neglected.
 * This always encloses the acceptance region accounting for electron mass without
 * compromising the acceptance efficiency.
 *
 * For a given \f$y\f$, the electron energy is then sampled
 * uniformly within the allowed Dalitz interval
 *
 * \f[ x_{\min} = 1-y+\frac{r^2}{1-y}, \qquad x_{\max} = 1+r^2, \f]
 *
 * reproducing the required joint distribution. All calculations are done in the rest
 * frame of the muon, with momentum vectors derived from E/p conservation. The final momentum
 * vectors are given an randomly sampled isotropic orientation in the muon rest frame
 * before returning the products.
 *
 * References:
 * - PDG 2026 Review, Muon Decay Parameters (Sec. 57)
 *   https://pdg.lbl.gov/2026/reviews/rpp2026-rev-muon-decay-params.pdf
 *   - \note The Geant4 implementation uses a different definition for reduced energies
 *     (\f$x\f$ and \f$y\f$) than the PDG. This gives slightly cleaner kinematic boundaries
 *     in code, and only leads to a difference in overall normalization factors.
 * - PDG 2026 Review, Kinematics, three-body decays and Dalitz form (Sec. 49.4.3)
 *   https://pdg.lbl.gov/2026/reviews/rpp2026-rev-kinematics.pdf
 *
 * \returns Heap-allocated \ref G4DecayProducts containing sampled decay products. Caller must take
 * ownership.
 */
G4DecayProducts* G4MuonDecayChannel::DecayIt(G4double)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) G4cout << "G4MuonDecayChannel::DecayIt ";
#endif

  CheckAndFillParent();
  CheckAndFillDaughters();

  // parent mass
  const G4double parentmass = G4MT_parent->GetPDGMass();
  const G4int N_DAUGHTER = 3;

  // daughters'mass
  G4double daughtermass[N_DAUGHTER];
  for (G4int index = 0; index < N_DAUGHTER; ++index)
  {
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
  }

  // create parent G4DynamicParticle at rest
  auto parentparticle = new G4DynamicParticle(G4MT_parent, G4ThreeVector(), 0.0);
  auto products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double x;  // Reduced electron energy.
  G4double y;  // Reduced electron-neutrino energy.

  const G4double Emax = 0.5 * parentmass;
  const G4double r = daughtermass[0] / parentmass;
  const G4double a = 1.0 - r * r;
  const G4double envelope_factor = 0.25 * 27.0;

  // Sample electron anti-neutrino and electron reduced energies
  while (true)
  {
    y = a * G4UniformRand();
    const G4double mass_correction_factor = a - y;
    const G4double f =
      envelope_factor * mass_correction_factor * mass_correction_factor * y * y / (1.0 - y);
    if (G4UniformRand() <= f) break;
  }

  const G4double x_min = 1.0 - y + r * r / (1.0 - y);
  const G4double x_max = 1.0 + r * r;
  const G4double dx = x_max - x_min;
  // Clamped to avoid negative dx at limits from precision
  x = (dx <= 0.0) ? x_max : x_min + dx * G4UniformRand();

  // Calculate physical energies and momenta
  const G4double Ee = Emax * x;
  const G4double Pe = std::sqrt(std::max(0.0, Ee * Ee - daughtermass[0] * daughtermass[0]));
  const G4double Enue = Emax * y;
  const G4double Enumu = parentmass - Ee - Enue;

  // Sample momentum vectors, handling singularities when one daughter particle is at rest.
  constexpr G4double kMomentumTolerance = 1e-9 * MeV;
  const G4ThreeVector random_direction = G4RandomDirection();
  G4ThreeVector direction0, direction1, direction2;

  if (Pe < kMomentumTolerance)
  {
    direction1 = Enue * random_direction;
    direction2 = -1.0 * Enumu * random_direction;
  }
  else if (Enue < kMomentumTolerance)
  {
    // Electron antineutrino at rest
    direction0 = Pe * random_direction;
    direction2 = -1.0 * Enumu * random_direction;
  }
  else if (Enumu < kMomentumTolerance)
  {
    // Muon neutrino at rest
    direction0 = Pe * random_direction;
    direction1 = -1.0 * Enue * random_direction;
  }
  else
  {
    // Kahan formula used to prevent numerical issues close to back-to-back
    // thresholds. Tangent-half-angle substitution to get cos/sin directly.
    // Note that Kahan gives the interior angle between Pe and Enue, so conversion
    // to _opening_ angle needed (phi = pi - theta, cos_theta = -cos_phi, sin_theta = sin_phi)
    const G4double tan_half_angle = kahan_interior_angle(Pe, Enue, Enumu);
    const G4double t_sq = tan_half_angle * tan_half_angle;
    const G4double denom = 1.0 / (1.0 + t_sq);

    const G4double cos_theta = (t_sq - 1.0) * denom;
    const G4double sin_theta = 2.0 * tan_half_angle * denom;

    // Place electron momentum along new z(polar)-axis
    direction0 = Pe * random_direction;

    // Construct electron-neutrino direction as if electron momentum along z-axis
    G4double rphi = twopi * G4UniformRand() * rad;
    direction1.set(Enue * sin_theta * std::cos(rphi), Enue * sin_theta * std::sin(rphi),
                   Enue * cos_theta);
    // Rotate to actual frame as defined by electron direction
    direction1.rotateUz(random_direction);

    direction2 = -1.0 * (direction0 + direction1);
  }

  // electron 0
  auto daughterparticle = new G4DynamicParticle(G4MT_daughters[0], direction0);
  products->PushProducts(daughterparticle);

  // electronic neutrino  1
  auto daughterparticle1 = new G4DynamicParticle(G4MT_daughters[1], direction1);
  products->PushProducts(daughterparticle1);

  // muonnic neutrino 2
  auto daughterparticle2 = new G4DynamicParticle(G4MT_daughters[2], direction2);
  products->PushProducts(daughterparticle2);

// output message
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1)
  {
    G4cout << "G4MuonDecayChannel::DecayIt()";
    G4cout << " create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }
#endif
  return products;
}
