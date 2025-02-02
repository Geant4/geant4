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

#ifndef G4FERMIDATATYPES_HH
#define G4FERMIDATATYPES_HH

#include <G4LorentzVector.hh>
#include <G4Vector3D.hh>
#include <G4String.hh>

#include <memory>
#include <string>

namespace fbu
{
using G4FermiInt = G4int;

using G4FermiUInt = uint32_t;

using G4FermiFloat = G4double;

using G4FermiLorentzVector = G4LorentzVector;

using G4FermiVector3 = G4Vector3D;

using G4FermiParticleMomentum = G4FermiVector3;

using G4FermiStr = G4String;

static constexpr G4FermiInt MAX_Z = 9;
static constexpr G4FermiInt MAX_A = 17;

template<typename Key, typename Value>
class G4FermiVCache
{
  public:
    virtual std::shared_ptr<Value> Insert(const Key& key, Value&& value) = 0;
    virtual std::shared_ptr<Value> Get(const Key& key) = 0;
    virtual ~G4FermiVCache() = default;
};

class G4FermiAtomicMass
{
  public:
    using G4FermiValueType = G4FermiUInt;

    G4FermiAtomicMass() = default;

    explicit constexpr G4FermiAtomicMass(G4FermiValueType mass) : mass_(mass) {}

    G4FermiAtomicMass(const G4FermiAtomicMass& other) = default;

    G4FermiAtomicMass(G4FermiAtomicMass&& other) = default;

    G4FermiAtomicMass& operator=(const G4FermiAtomicMass& other) = default;

    G4FermiAtomicMass& operator=(G4FermiAtomicMass&& other) = default;

    constexpr operator G4FermiUInt() const { return mass_; }

    constexpr operator G4FermiInt() const { return mass_; }

    constexpr operator G4FermiFloat() const { return mass_; }

    bool operator<(const G4FermiAtomicMass& other) const { return mass_ < other.mass_; }

    bool operator>(const G4FermiAtomicMass& other) const { return mass_ > other.mass_; }

    bool operator<=(const G4FermiAtomicMass& other) const { return mass_ <= other.mass_; }

    bool operator>=(const G4FermiAtomicMass& other) const { return mass_ >= other.mass_; }

    bool operator==(const G4FermiAtomicMass& other) const { return mass_ == other.mass_; }

    bool operator!=(const G4FermiAtomicMass& other) const { return mass_ != other.mass_; }

  private:
    G4FermiValueType mass_;
};

class G4FermiChargeNumber
{
  public:
    using G4FermiValueType = G4FermiUInt;

    G4FermiChargeNumber() = default;

    explicit constexpr G4FermiChargeNumber(G4FermiValueType charge) : charge_(charge) {}

    G4FermiChargeNumber(const G4FermiChargeNumber& other) = default;

    G4FermiChargeNumber(G4FermiChargeNumber&& other) = default;

    G4FermiChargeNumber& operator=(const G4FermiChargeNumber& other) = default;

    G4FermiChargeNumber& operator=(G4FermiChargeNumber&& other) = default;

    constexpr operator G4FermiUInt() const { return charge_; }

    constexpr operator G4FermiInt() const { return charge_; }

    constexpr operator G4FermiFloat() const { return charge_; }

    bool operator<(const G4FermiChargeNumber& other) const { return charge_ < other.charge_; }

    bool operator>(const G4FermiChargeNumber& other) const { return charge_ > other.charge_; }

    bool operator<=(const G4FermiChargeNumber& other) const { return charge_ <= other.charge_; }

    bool operator>=(const G4FermiChargeNumber& other) const { return charge_ >= other.charge_; }

    bool operator==(const G4FermiChargeNumber& other) const { return charge_ == other.charge_; }

    bool operator!=(const G4FermiChargeNumber& other) const { return charge_ != other.charge_; }

  private:
    G4FermiValueType charge_;
};

struct G4FermiNucleiData
{
    G4FermiAtomicMass atomicMass;
    G4FermiChargeNumber chargeNumber;

    bool operator<(const G4FermiNucleiData& other) const
    {
      return atomicMass < other.atomicMass
             || (atomicMass == other.atomicMass && chargeNumber < other.chargeNumber);
    }

    bool operator==(const G4FermiNucleiData& other) const
    {
      return atomicMass == other.atomicMass && chargeNumber == other.chargeNumber;
    }

    bool operator!=(const G4FermiNucleiData& other) const
    {
      return atomicMass != other.atomicMass || chargeNumber != other.chargeNumber;
    }
};

}  // namespace fbu

namespace std
{
template<>
struct hash<::fbu::G4FermiNucleiData>
{
    size_t operator()(const ::fbu::G4FermiNucleiData& key) const
    {
      auto mass = ::fbu::G4FermiInt(key.atomicMass);
      auto charge = ::fbu::G4FermiInt(key.chargeNumber);
      return (mass * (mass + 1)) / 2 + charge;
    }
};

string to_string(::fbu::G4FermiAtomicMass mass);
string to_string(::fbu::G4FermiChargeNumber charge);

std::ostream& operator<<(std::ostream& out, const ::fbu::G4FermiAtomicMass& mass);
std::istream& operator>>(std::istream& in, ::fbu::G4FermiAtomicMass& mass);

std::ostream& operator<<(std::ostream& out, const ::fbu::G4FermiChargeNumber& charge);
std::istream& operator>>(std::istream& in, ::fbu::G4FermiChargeNumber& charge);
}  // namespace std

constexpr fbu::G4FermiAtomicMass operator""_m(unsigned long long mass)
{
  return fbu::G4FermiAtomicMass(mass);
}

constexpr fbu::G4FermiChargeNumber operator""_c(unsigned long long charge)
{
  return fbu::G4FermiChargeNumber(charge);
}

#endif  // G4FERMIDATATYPES_HH
