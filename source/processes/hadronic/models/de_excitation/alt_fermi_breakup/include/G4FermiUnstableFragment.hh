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

#ifndef G4FERMIUNSTABLEFRAGMENT_HH
#define G4FERMIUNSTABLEFRAGMENT_HH

#include "G4FermiPossibleFragment.hh"

namespace fbu
{

class G4FermiUnstableFragment : public G4FermiPossibleFragment
{
  public:
    G4FermiUnstableFragment(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                            G4FermiInt polarization, G4FermiFloat excitationEnergy,
                            std::vector<G4FermiNucleiData>&& decayData);

    void AppendDecayFragments(const G4FermiLorentzVector& momentum,
                              std::vector<G4FermiParticle>& particles) const override;

  private:
    std::vector<G4FermiNucleiData> decayData_;
    std::vector<G4FermiFloat> masses_;
};

#define FERMI_ADD_UNSTABLE_FRAGMENT(NAME, FRAGMENTS)                                                      \
  inline G4FermiUnstableFragment NAME(G4FermiAtomicMass atomicMass,                              \
                                      G4FermiChargeNumber chargeNumber, G4FermiInt polarization, \
                                      G4FermiFloat excitationEnergy)                             \
  {                                                                                              \
    return G4FermiUnstableFragment(atomicMass, chargeNumber, polarization, excitationEnergy,     \
                                   FRAGMENTS);                                                   \
  }

// He5 ----> alpha + neutron
FERMI_ADD_UNSTABLE_FRAGMENT(He5Fragment, std::vector<G4FermiNucleiData>({
                                  G4FermiNucleiData{4_m, 2_c},
                                  G4FermiNucleiData{1_m, 0_c},
                                }));

// B9 ----> alpha + alpha + proton
FERMI_ADD_UNSTABLE_FRAGMENT(B9Fragment, std::vector<G4FermiNucleiData>({
                                 G4FermiNucleiData{4_m, 2_c},
                                 G4FermiNucleiData{4_m, 2_c},
                                 G4FermiNucleiData{1_m, 1_c},
                               }));

// Be8 ----> alpha + alpha
FERMI_ADD_UNSTABLE_FRAGMENT(Be8Fragment, std::vector<G4FermiNucleiData>({
                                  G4FermiNucleiData{4_m, 2_c},
                                  G4FermiNucleiData{4_m, 2_c},
                                }));

// Li5 ----> alpha + proton
FERMI_ADD_UNSTABLE_FRAGMENT(Li5Fragment, std::vector<G4FermiNucleiData>({
                                  G4FermiNucleiData{4_m, 2_c},
                                  G4FermiNucleiData{1_m, 1_c},
                                }));

#undef FERMI_ADD_UNSTABLE_FRAGMENT

}  // namespace fbu

#endif  // G4FERMIUNSTABLEFRAGMENT_HH
