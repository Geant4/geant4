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

#ifndef G4FERMIDATATYPES_HH
#define G4FERMIDATATYPES_HH

#include <G4LorentzVector.hh>
#include <G4PhysicalConstants.hh>
#include <G4String.hh>
#include <G4Vector3D.hh>
#include <Randomize.hh>
#include <globals.hh>

#include <G4qss_misc.hh>

static constexpr G4int MAX_Z = 9;
static constexpr G4int MAX_A = 17;

G4Vector3D SampleIsotropicVector(G4double magnitude);

class G4FermiAtomicMass
{
  public:
    using ValueType = std::uint32_t;

    G4FermiAtomicMass() = default;

    explicit constexpr G4FermiAtomicMass(ValueType mass) : mass_(mass) {}

    G4FermiAtomicMass(const G4FermiAtomicMass& other) = default;

    G4FermiAtomicMass(G4FermiAtomicMass&& other) = default;

    G4FermiAtomicMass& operator=(const G4FermiAtomicMass& other) = default;

    G4FermiAtomicMass& operator=(G4FermiAtomicMass&& other) = default;

    constexpr operator std::uint32_t() const { return mass_; }

    constexpr operator G4int() const { return mass_; }

    constexpr operator G4double() const { return mass_; }

    G4bool operator<(const G4FermiAtomicMass& other) const { return mass_ < other.mass_; }

    G4bool operator>(const G4FermiAtomicMass& other) const { return mass_ > other.mass_; }

    G4bool operator<=(const G4FermiAtomicMass& other) const { return mass_ <= other.mass_; }

    G4bool operator>=(const G4FermiAtomicMass& other) const { return mass_ >= other.mass_; }

    G4bool operator==(const G4FermiAtomicMass& other) const { return mass_ == other.mass_; }

    G4bool operator!=(const G4FermiAtomicMass& other) const { return mass_ != other.mass_; }

  private:
    ValueType mass_;
};

class G4FermiChargeNumber
{
  public:
    using ValueType = std::uint32_t;

    G4FermiChargeNumber() = default;

    explicit constexpr G4FermiChargeNumber(ValueType charge) : charge_(charge) {}

    G4FermiChargeNumber(const G4FermiChargeNumber& other) = default;

    G4FermiChargeNumber(G4FermiChargeNumber&& other) = default;

    G4FermiChargeNumber& operator=(const G4FermiChargeNumber& other) = default;

    G4FermiChargeNumber& operator=(G4FermiChargeNumber&& other) = default;

    constexpr operator std::uint32_t() const { return charge_; }

    constexpr operator G4int() const { return charge_; }

    constexpr operator G4double() const { return charge_; }

    G4bool operator<(const G4FermiChargeNumber& other) const { return charge_ < other.charge_; }

    G4bool operator>(const G4FermiChargeNumber& other) const { return charge_ > other.charge_; }

    G4bool operator<=(const G4FermiChargeNumber& other) const { return charge_ <= other.charge_; }

    G4bool operator>=(const G4FermiChargeNumber& other) const { return charge_ >= other.charge_; }

    G4bool operator==(const G4FermiChargeNumber& other) const { return charge_ == other.charge_; }

    G4bool operator!=(const G4FermiChargeNumber& other) const { return charge_ != other.charge_; }

  private:
    ValueType charge_;
};

struct G4FermiNucleiData
{
    G4FermiAtomicMass atomicMass;
    G4FermiChargeNumber chargeNumber;

    G4bool operator<(const G4FermiNucleiData& other) const
    {
      return atomicMass < other.atomicMass
             || (atomicMass == other.atomicMass && chargeNumber < other.chargeNumber);
    }

    G4bool operator==(const G4FermiNucleiData& other) const
    {
      return atomicMass == other.atomicMass && chargeNumber == other.chargeNumber;
    }

    G4bool operator!=(const G4FermiNucleiData& other) const
    {
      return atomicMass != other.atomicMass || chargeNumber != other.chargeNumber;
    }
};

namespace std
{
template<>
struct hash<G4FermiNucleiData>
{
    std::size_t operator()(const G4FermiNucleiData& key) const
    {
      auto mass = G4int(key.atomicMass);
      auto charge = G4int(key.chargeNumber);
      return (mass * (mass + 1)) / 2 + charge;
    }
};

string to_string(G4FermiAtomicMass mass);
string to_string(G4FermiChargeNumber charge);

std::ostream& operator<<(std::ostream& out, const G4FermiAtomicMass& mass);
std::istream& operator>>(std::istream& in, G4FermiAtomicMass& mass);

std::ostream& operator<<(std::ostream& out, const G4FermiChargeNumber& charge);
std::istream& operator>>(std::istream& in, G4FermiChargeNumber& charge);
}  // namespace std

constexpr G4FermiAtomicMass operator""_m(unsigned long long mass)
{
  return G4FermiAtomicMass(mass);
}

constexpr G4FermiChargeNumber operator""_c(unsigned long long charge)
{
  return G4FermiChargeNumber(charge);
}

#define FERMI_ASSERT_MSG(COND, MSG)                                                             \
  if (unlikely(!(COND))) {                                                                      \
    std::ostringstream sstream;                                                                 \
    sstream << "assertion failed: \"" << #COND << '\"' << " at " << __FILE__ << ':' << __LINE__ \
            << '\n'                                                                             \
            << MSG;                                                                             \
    throw std::runtime_error(sstream.str());                                                    \
  }

#endif  // G4FERMIDATATYPES_HH
