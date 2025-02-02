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

#include "G4FermiParticle.hh"

#include "G4FermiLogger.hh"
#include "G4FermiNucleiProperties.hh"

#include <CLHEP/Units/PhysicalConstants.h>

#include <iomanip>

using namespace fbu;

G4FermiParticle::G4FermiParticle(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                                 const G4FermiLorentzVector& momentum)
  : atomicMass_(atomicMass),
    chargeNumber_(chargeNumber),
    momentum_(momentum),
    groundStateMass_(G4FermiNucleiProperties()->GetNuclearMass(atomicMass_, chargeNumber_))
{
  FERMI_ASSERT_MSG(G4FermiUInt(atomicMass_) >= G4FermiUInt(chargeNumber),
                   "imposible particle: A = " << atomicMass_ << ", Z = " << chargeNumber);

  CalculateExcitationEnergy();
}

G4FermiNucleiData G4FermiParticle::GetNucleiData() const
{
  return {
    atomicMass_,
    chargeNumber_,
  };
}

G4FermiAtomicMass G4FermiParticle::GetAtomicMass() const
{
  return atomicMass_;
}

G4FermiChargeNumber G4FermiParticle::GetChargeNumber() const
{
  return chargeNumber_;
}

const G4FermiLorentzVector& G4FermiParticle::GetMomentum() const
{
  return momentum_;
}

G4FermiFloat G4FermiParticle::GetExcitationEnergy() const
{
  return excitationEnergy_;
}

G4FermiFloat G4FermiParticle::GetGroundStateMass() const
{
  return groundStateMass_;
}

bool G4FermiParticle::IsStable() const
{
  return excitationEnergy_ <= 0.;
}

void G4FermiParticle::CalculateExcitationEnergy()
{
  excitationEnergy_ = momentum_.mag() - groundStateMass_;
  if (excitationEnergy_ < 0.) {
    if (excitationEnergy_ < -10 * CLHEP::eV) {
      FERMI_LOG_WARN("Excitation Energy is too negative: " << excitationEnergy_ / CLHEP::MeV
                                                           << " MeV");
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
