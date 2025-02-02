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

#ifndef G4FERMIPOSSIBLEFRAGMENT_HH
#define G4FERMIPOSSIBLEFRAGMENT_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiParticle.hh"

#include <ostream>
#include <vector>

namespace fbu
{
class G4FermiPossibleFragment;

using G4FermiPossibleFragmentVector = std::vector<const G4FermiPossibleFragment*>;

class G4FermiPossibleFragment
{
  public:
    G4FermiPossibleFragment(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                    G4FermiInt polarization, G4FermiFloat excitationEnergy);

    G4FermiPossibleFragment(const G4FermiPossibleFragment&) = delete;

    G4FermiPossibleFragment& operator=(const G4FermiPossibleFragment&) = delete;

    std::vector<G4FermiParticle> GetDecayFragments(const G4FermiLorentzVector& momentum) const;

    virtual void AppendDecayFragments(const G4FermiLorentzVector& momentum,
                                      std::vector<G4FermiParticle>& particles) const = 0;

    G4FermiAtomicMass GetAtomicMass() const;

    G4FermiChargeNumber GetChargeNumber() const;

    G4FermiInt GetPolarization() const;

    G4FermiFloat GetExcitationEnergy() const;

    G4FermiFloat GetMass() const;

    G4FermiFloat GetTotalEnergy() const;

    virtual ~G4FermiPossibleFragment() = default;

  protected:
    G4FermiAtomicMass atomicMass_;  // A
    G4FermiChargeNumber chargeNumber_;  // Z
    G4FermiInt polarization_;
    G4FermiFloat excitationEnergy_;
};

}  // namespace fbu

namespace std
{
ostream& operator<<(ostream&, const ::fbu::G4FermiPossibleFragment&);
}  // namespace std
#endif  // G4FERMIPOSSIBLEFRAGMENT_HH
