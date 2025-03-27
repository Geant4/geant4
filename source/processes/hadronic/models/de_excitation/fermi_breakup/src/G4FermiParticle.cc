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

#include "G4FermiParticle.hh"

#include "G4FermiDataTypes.hh"
#include "G4FermiNucleiProperties.hh"

#include <G4PhysicalConstants.hh>

#include <iomanip>

G4FermiParticle::G4FermiParticle(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                                 const G4LorentzVector& momentum)
  : atomicMass_(atomicMass), chargeNumber_(chargeNumber), momentum_(momentum)
{
  FERMI_ASSERT_MSG(static_cast<std::uint32_t>(atomicMass_)
                     >= static_cast<std::uint32_t>(chargeNumber),
                   "imposible particle: A = " << atomicMass_ << ", Z = " << chargeNumber);

  RecalculateExcitationEnergy();
}

G4FermiAtomicMass G4FermiParticle::GetAtomicMass() const
{
  return atomicMass_;
}

G4FermiChargeNumber G4FermiParticle::GetChargeNumber() const
{
  return chargeNumber_;
}

const G4LorentzVector& G4FermiParticle::GetMomentum() const
{
  return momentum_;
}

G4double G4FermiParticle::GetExcitationEnergy() const
{
  return excitationEnergy_;
}

G4bool G4FermiParticle::IsStable() const
{
  return excitationEnergy_ <= 0.;
}

void G4FermiParticle::RecalculateExcitationEnergy()
{
  excitationEnergy_ =
    momentum_.mag() - G4FermiNucleiProperties::GetNuclearMass(atomicMass_, chargeNumber_);
  if (excitationEnergy_ < 0.) {
    if (excitationEnergy_ < -10.0 * CLHEP::eV) {
      G4ExceptionDescription ed;
      ed << "Excitation energy is too negative: " << excitationEnergy_ / CLHEP::MeV << " MeV";
      G4Exception("G4FermiParticle::RecalculateExcitationEnergy()", "Fermi001", JustWarning, ed);
    }
    excitationEnergy_ = 0.;
  }
}

std::ostream& std::operator<<(std::ostream& out, const G4FermiParticle& particle)
{
  const auto oldFlags = out.flags();
  const auto oldUserPrecision = out.precision();

  out.setf(std::ios::floatfield);
  out << "FermiParticle: { A = " << particle.GetAtomicMass()
      << ", Z = " << particle.GetChargeNumber();

  out.setf(std::ios::scientific, std::ios::floatfield);
  out << std::setprecision(3) << ", U = " << particle.GetExcitationEnergy() / CLHEP::MeV << " MeV"
      << ", IsGroundState = " << (particle.IsStable() ? "yes" : "no") << ", P = ("
      << particle.GetMomentum().x() / CLHEP::MeV << ", " << particle.GetMomentum().y() / CLHEP::MeV
      << ", " << particle.GetMomentum().z() / CLHEP::MeV
      << ") MeV, E = " << particle.GetMomentum().t() / CLHEP::MeV << " MeV}"
      << " }";

  out.setf(oldFlags, std::ios::floatfield);
  out.precision(oldUserPrecision);

  return out;
}
