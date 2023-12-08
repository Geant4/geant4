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
// G4VDecayChannel
//
// Class description:
//
// Abstract class to describe decay kinematics

// Author: H.Kurashige, 27 July 1996
// --------------------------------------------------------------------
#ifndef G4VDecayChannel_hh
#define G4VDecayChannel_hh 1

#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <cmath>

class G4ParticleDefinition;
class G4DecayProducts;
class G4ParticleTable;

class G4VDecayChannel
{
  public:
    // Constructors
    G4VDecayChannel(const G4String& aName, G4int Verbose = 1);
    G4VDecayChannel(const G4String& aName, const G4String& theParentName, G4double theBR,
                    G4int theNumberOfDaughters, const G4String& theDaughterName1,
                    const G4String& theDaughterName2 = "", const G4String& theDaughterName3 = "",
                    const G4String& theDaughterName4 = "", const G4String& theDaughterName5 = "");

    // Destructor
    virtual ~G4VDecayChannel();

    // Equality operators
    G4bool operator==(const G4VDecayChannel& r) const { return (this == &r); }
    G4bool operator!=(const G4VDecayChannel& r) const { return (this != &r); }

    // Less-than operator is defined for G4DecayTable
    inline G4bool operator<(const G4VDecayChannel& right) const;

    virtual G4DecayProducts* DecayIt(G4double parentMass = -1.0) = 0;

    // Get kinematics name
    inline const G4String& GetKinematicsName() const;

    // Get branching ratio
    inline G4double GetBR() const;

    // Get number of daughter particles
    inline G4int GetNumberOfDaughters() const;

    // Get the pointer to the parent particle
    inline G4ParticleDefinition* GetParent();

    // Get the pointer to a daughter particle
    inline G4ParticleDefinition* GetDaughter(G4int anIndex);

    // Get the angular momentum of the decay
    G4int GetAngularMomentum();

    // Get the name of the parent particle
    inline const G4String& GetParentName() const;

    // Get the name of a daughter particle
    inline const G4String& GetDaughterName(G4int anIndex) const;

    // Get mass of parent
    inline G4double GetParentMass() const;

    // Get mass of daughter
    inline G4double GetDaughterMass(G4int anIndex) const;

    // Set the parent particle (by name or by pointer)
    void SetParent(const G4ParticleDefinition* particle_type);
    inline void SetParent(const G4String& particle_name);

    // Set branching ratio
    void SetBR(G4double value);

    // Set number of daughter particles
    void SetNumberOfDaughters(G4int value);

    // Set a daughter particle (by name or by pointer)
    void SetDaughter(G4int anIndex, const G4ParticleDefinition* particle_type);
    void SetDaughter(G4int anIndex, const G4String& particle_name);

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;
    void DumpInfo();

    inline G4double GetRangeMass() const;
    inline void SetRangeMass(G4double val);
    virtual G4bool IsOKWithParentMass(G4double parentMass);

    void SetPolarization(const G4ThreeVector&);
    inline const G4ThreeVector& GetPolarization() const;

  protected:
    // Default constructor
    G4VDecayChannel();

    // Copy constructor and assignment operator
    G4VDecayChannel(const G4VDecayChannel&);
    G4VDecayChannel& operator=(const G4VDecayChannel&);

    // Clear daughters array
    void ClearDaughtersName();

    inline void CheckAndFillDaughters();
    inline void CheckAndFillParent();

    G4double DynamicalMass(G4double massPDG, G4double width, G4double maxDev = 1.0) const;

  protected:
    // Kinematics name
    G4String kinematics_name = "";
    // Branching ratio  [0.0 - 1.0]
    G4double rbranch = 0.0;
    // Parent particle
    G4String* parent_name = nullptr;
    // Daughter particles
    G4String** daughters_name = nullptr;

    // Range of mass allowed in decay
    G4double rangeMass = 2.5;

    // Polarization of the parent particle
    G4ThreeVector parent_polarization;

    // Pointer to particle table
    G4ParticleTable* particletable = nullptr;

    static const G4String noName;

    G4ParticleDefinition* G4MT_parent = nullptr;
    G4ParticleDefinition** G4MT_daughters = nullptr;
    G4double G4MT_parent_mass = 0.0;
    G4double* G4MT_daughters_mass = nullptr;
    G4double* G4MT_daughters_width = nullptr;
    G4Mutex daughtersMutex;
    G4Mutex parentMutex;

    // Number of daughters
    G4int numberOfDaughters = 0;

    // Control flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    G4int verboseLevel = 1;

  private:
    // Fill daughters array
    void FillDaughters();

    // Fill parent
    void FillParent();

    const G4String& GetNoName() const;
};

// ------------------------------------------------------------
// Inline methods
// ------------------------------------------------------------

inline G4bool G4VDecayChannel::operator<(const G4VDecayChannel& right) const
{
  return (this->rbranch < right.rbranch);
}

inline G4ParticleDefinition* G4VDecayChannel::GetDaughter(G4int anIndex)
{
  // pointers to daughter particles are filled, if they are not set yet
  CheckAndFillDaughters();

  // get the pointer to a daughter particle
  if ((anIndex >= 0) && (anIndex < numberOfDaughters)) {
    return G4MT_daughters[anIndex];
  }

  if (verboseLevel > 0)
    G4cout << "G4VDecayChannel::GetDaughter  index out of range " << anIndex << G4endl;
  return nullptr;
}

inline const G4String& G4VDecayChannel::GetDaughterName(G4int anIndex) const
{
  if ((anIndex >= 0) && (anIndex < numberOfDaughters)) {
    return *daughters_name[anIndex];
  }

  if (verboseLevel > 0) {
    G4cout << "G4VDecayChannel::GetDaughterName ";
    G4cout << "index out of range " << anIndex << G4endl;
  }
  return GetNoName();
}

inline G4double G4VDecayChannel::GetDaughterMass(G4int anIndex) const
{
  if ((anIndex >= 0) && (anIndex < numberOfDaughters)) {
    return G4MT_daughters_mass[anIndex];
  }

  if (verboseLevel > 0) {
    G4cout << "G4VDecayChannel::GetDaughterMass ";
    G4cout << "index out of range " << anIndex << G4endl;
  }
  return 0.0;
}

inline G4ParticleDefinition* G4VDecayChannel::GetParent()
{
  // the pointer to the parent particle is filled, if it is not set yet
  CheckAndFillParent();
  // get the pointer to the parent particle
  return G4MT_parent;
}

inline const G4String& G4VDecayChannel::GetParentName() const
{
  return *parent_name;
}

inline G4double G4VDecayChannel::GetParentMass() const
{
  return G4MT_parent_mass;
}

inline void G4VDecayChannel::SetParent(const G4String& particle_name)
{
  delete parent_name;
  parent_name = new G4String(particle_name);
  G4MT_parent = nullptr;
}

inline G4int G4VDecayChannel::GetNumberOfDaughters() const
{
  return numberOfDaughters;
}

inline const G4String& G4VDecayChannel::GetKinematicsName() const
{
  return kinematics_name;
}

inline G4double G4VDecayChannel::GetBR() const
{
  return rbranch;
}

inline void G4VDecayChannel::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline G4int G4VDecayChannel::GetVerboseLevel() const
{
  return verboseLevel;
}

inline G4double G4VDecayChannel::GetRangeMass() const
{
  return rangeMass;
}

inline void G4VDecayChannel::SetRangeMass(G4double val)
{
  if (val >= 0.) rangeMass = val;
}

inline void G4VDecayChannel::SetPolarization(const G4ThreeVector& polar)
{
  parent_polarization = polar;
}

inline const G4ThreeVector& G4VDecayChannel::GetPolarization() const
{
  return parent_polarization;
}

inline void G4VDecayChannel::CheckAndFillDaughters()
{
  G4AutoLock l(&daughtersMutex);
  if (G4MT_daughters == nullptr) {
    l.unlock();
    FillDaughters();
  }
}

inline void G4VDecayChannel::CheckAndFillParent()
{
  G4AutoLock l(&parentMutex);
  if (G4MT_parent == nullptr) {
    l.unlock();
    FillParent();
  }
}

#endif
