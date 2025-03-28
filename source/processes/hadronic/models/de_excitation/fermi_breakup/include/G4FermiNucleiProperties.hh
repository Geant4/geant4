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

#ifndef G4FERMINUCLEIPROPERTIES_HH
#define G4FERMINUCLEIPROPERTIES_HH

#include "G4FermiDataTypes.hh"

// Caches values from larger G4NucleiProperties(x5-10 speed boost)
class G4FermiNucleiProperties
{
  public:
    void Initialize() { *this = G4FermiNucleiProperties(); }

    template<typename DataSource>
    void Initialize(const DataSource& dataSource)
    {
      Initialize(dataSource.begin(), dataSource.end());
    }

    template<typename Iter>
    void Initialize(Iter begin, Iter end)
    {
      nucleiMasses_.clear();
      static_assert(
        std::is_same_v<typename Iter::value_type, std::pair<const G4FermiNucleiData, G4double>>,
        "invalid iterator");
      for (auto it = begin; it != end; ++it) {
        InsertNuclei(it->first.atomicMass, it->first.chargeNumber, it->second);
      }
    }

    static G4double GetNuclearMass(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
    {
      return Instance().GetNuclearMassImpl(atomicMass, chargeNumber);
    }

    static G4bool IsStable(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber)
    {
      return Instance().IsStableImpl(atomicMass, chargeNumber);
    }

    void InsertNuclei(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber, G4double mass,
                      G4bool isStable = true);

    static G4FermiNucleiProperties& Instance()
    {
      static G4FermiNucleiProperties properties;
      return properties;
    }

  private:
    G4FermiNucleiProperties();

    G4double GetNuclearMassImpl(G4FermiAtomicMass atomicMass,
                                G4FermiChargeNumber chargeNumber) const;

    G4bool IsStableImpl(G4FermiAtomicMass atomicMass, G4FermiChargeNumber chargeNumber) const;

    struct G4FermiMassData
    {
        G4double mass;

        G4bool isStable = false;  // is nuclei stable

        G4bool isCached = false;  // value has been inserted earlier
    };

    mutable std::vector<G4FermiMassData> nucleiMasses_;
};

#endif  // G4FERMINUCLEIPROPERTIES_HH
