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

#ifndef G4FERMIPARTICLE_HH
#define G4FERMIPARTICLE_HH

#include "G4FermiDataTypes.hh"

namespace fbu
{

class G4FermiParticle
{
  public:
    G4FermiParticle() = delete;

    G4FermiParticle(const G4FermiParticle&) = default;
    G4FermiParticle(G4FermiParticle&&) = default;

    G4FermiParticle& operator=(const G4FermiParticle&) = default;
    G4FermiParticle& operator=(G4FermiParticle&&) = default;

    G4FermiParticle(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                    const G4FermiLorentzVector& momentum);

    G4FermiNucleiData GetNucleiData() const;

    G4FermiAtomicMass GetAtomicMass() const;

    G4FermiChargeNumber GetChargeNumber() const;

    const G4FermiLorentzVector& GetMomentum() const;

    G4FermiFloat GetExcitationEnergy() const;

    G4FermiFloat GetGroundStateMass() const;

    bool IsStable() const;

  private:
    void CalculateExcitationEnergy();

    G4FermiAtomicMass atomicMass_;
    G4FermiChargeNumber chargeNumber_;
    G4FermiLorentzVector momentum_;

    G4FermiFloat groundStateMass_ = 0;
    G4FermiFloat excitationEnergy_ = 0;
};

}  // namespace fbu

namespace std
{
ostream& operator<<(ostream&, const ::fbu::G4FermiParticle&);
}  // namespace std

#endif  // G4FERMIPARTICLE_HH
