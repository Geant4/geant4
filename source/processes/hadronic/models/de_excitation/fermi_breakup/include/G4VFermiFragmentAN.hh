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
// G4FermiBreakUpAN alternative FermiBreakUp model
// by A. Novikov (January 2025)
//

#ifndef G4VFERMIFRAGMENTAN_HH
#define G4VFERMIFRAGMENTAN_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiParticle.hh"

#include "globals.hh"
#include <vector>

class G4VFermiFragmentAN;

using G4FermiFragmentVector = std::vector<const G4VFermiFragmentAN*>;

class G4VFermiFragmentAN
{
  public:
    G4VFermiFragmentAN(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                       G4int polarization, G4double excitationEnergy);

    G4VFermiFragmentAN(const G4VFermiFragmentAN&) = delete;

    G4VFermiFragmentAN& operator=(const G4VFermiFragmentAN&) = delete;

    ~G4VFermiFragmentAN() = default;
  
    void Initialize();

    std::vector<G4FermiParticle> GetDecayFragments(const G4LorentzVector& momentum) const;

    virtual void AppendDecayFragments(const G4LorentzVector& momentum,
                                      std::vector<G4FermiParticle>& particles) const = 0;

    G4FermiAtomicMass GetAtomicMass() const;

    G4FermiChargeNumber GetChargeNumber() const;

    G4int GetPolarization() const;

    G4double GetExcitationEnergy() const;

    G4double GetMass() const;

    G4double GetTotalEnergy() const;

  protected:
    virtual void DoInitialize() = 0;

    G4FermiAtomicMass atomicMass_;  // A
    G4FermiChargeNumber chargeNumber_;  // Z
    G4int polarization_;

    G4double groudStateMass_;
    G4double excitationEnergy_;
};

namespace std
{
ostream& operator<<(ostream&, const G4VFermiFragmentAN&);
}  // namespace std

#endif  // G4VFERMIFRAGMENTAN_HH
