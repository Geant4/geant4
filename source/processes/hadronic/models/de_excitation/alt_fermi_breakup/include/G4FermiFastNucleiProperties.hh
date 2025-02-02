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

#ifndef G4FERMIFASTNUCLEIPROPERTIES_HH
#define G4FERMIFASTNUCLEIPROPERTIES_HH

#include "G4FermiDataTypes.hh"
#include "G4FermiVNucleiProperties.hh"

#include <vector>

namespace fbu
{
class G4FermiFastNucleiProperties : public G4FermiVNucleiProperties
{
  public:
    G4FermiFastNucleiProperties();

    template<typename DataSource>
    G4FermiFastNucleiProperties(const DataSource& dataSource);

    template<typename Iter>
    G4FermiFastNucleiProperties(Iter begin, Iter end);

    [[nodiscard]] G4FermiFloat GetNuclearMass(G4FermiAtomicMass atomicMass,
                                              G4FermiChargeNumber chargeNumber) const override;

    [[nodiscard]] bool IsStable(G4FermiAtomicMass atomicMass,
                                G4FermiChargeNumber chargeNumber) const override;

    void AddStableNuclei(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber,
                         G4FermiFloat mass);

    void AddStableNuclei(G4FermiNucleiData nucleiData, G4FermiFloat mass);

  private:
    struct G4FermiMassData
    {
        G4FermiFloat mass;

        bool isStable = false;  // value was added via AddMass, it is considered stable

        bool isCached = false;  // value has been calculated earlier
    };

    mutable std::vector<G4FermiMassData> nucleiMasses_;
};

template<typename DataSource>
G4FermiFastNucleiProperties::G4FermiFastNucleiProperties(const DataSource& dataSource)
  : G4FermiFastNucleiProperties(dataSource.begin(), dataSource.end())
{}

template<typename Iter>
G4FermiFastNucleiProperties::G4FermiFastNucleiProperties(Iter begin, Iter end)
{
  static_assert(
    std::is_same_v<typename Iter::value_type, std::pair<const G4FermiNucleiData, G4FermiFloat>>,
    "invalid iterator");
  for (auto it = begin; it != end; ++it) {
    AddStableNuclei(it->first, it->second);
  }
}

}  // namespace fbu

#endif  // G4FERMIFASTNUCLEIPROPERTIES_HH
