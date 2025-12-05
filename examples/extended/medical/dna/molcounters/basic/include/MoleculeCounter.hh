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
/// \file MoleculeCounter.hh
/// \brief Definition of the MoleculeCounter class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#ifndef MoleculeCounter_hh
#define MoleculeCounter_hh 1

#include "G4MolecularConfiguration.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4VUserMoleculeCounter.hh"

#include <map>
#include <set>
#include <vector>

class G4Navigator;

class G4StepPoint;

class G4Material;

class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct MoleculeCounterIndex : G4VMoleculeCounter::G4VMoleculeCounterIndex
{
    const G4MolecularConfiguration* Molecule;
    const G4VPhysicalVolume* Volume;
    std::vector<G4int> CopyNumbers;

    MoleculeCounterIndex() : Molecule(nullptr), Volume(nullptr), CopyNumbers() {}

    MoleculeCounterIndex(const MoleculeCounterIndex&) = default;

    MoleculeCounterIndex(MoleculeCounterIndex&&) = default;

    MoleculeCounterIndex(const G4MolecularConfiguration* molecule, const G4VTouchable* touchable)
      : Molecule(molecule), CopyNumbers()
    {
      if (touchable == nullptr) {
        Volume = nullptr;
      }
      else {
        Volume = touchable->GetVolume();
        CopyNumbers.reserve(touchable->GetHistoryDepth());
        for (auto i = 0; i < touchable->GetHistoryDepth(); ++i)
          CopyNumbers.push_back(touchable->GetCopyNumber(i));
      }
    }

    MoleculeCounterIndex(const G4MolecularConfiguration* molecule,
                         const G4VPhysicalVolume* volume,  // should be removed ?
                         const std::vector<G4int>& copyNumbers)
      : Molecule(molecule), Volume(volume), CopyNumbers(copyNumbers)
    {}

    ~MoleculeCounterIndex() override = default;

    MoleculeCounterIndex& operator=(const MoleculeCounterIndex&) = default;

    MoleculeCounterIndex& operator=(MoleculeCounterIndex&&) = default;

    G4bool operator<(G4VMoleculeCounterIndex const& other) const override
    {
      return *this < dynamic_cast<MoleculeCounterIndex const&>(other);
    }

    G4bool operator<(MoleculeCounterIndex const& other) const
    {
      if (std::less<>{}(Molecule, other.Molecule)) return true;
      if (std::less<>{}(other.Molecule, Molecule)) return false;

      if (std::less<>{}(Volume, other.Volume)) return true;
      if (std::less<>{}(other.Volume, Volume)) return false;

      return CopyNumbers < other.CopyNumbers;
    }

    G4bool operator==(G4VMoleculeCounterIndex const& other) const override
    {
      return *this == dynamic_cast<MoleculeCounterIndex const&>(other);
    }

    G4bool operator==(const MoleculeCounterIndex& other) const
    {
      return std::tie(Molecule, Volume) == std::tie(other.Molecule, other.Volume)
             && CopyNumbers == other.CopyNumbers;
    }

    G4String GetInfo() const override
    {
      G4String info = "(";
      const G4String& null = "null";
      info += Molecule == nullptr ? null : Molecule->GetName();
      info += ", ";
      info += Volume == nullptr ? null : Volume->GetName();
      info += ")";
      return info;
    }

    const G4MolecularConfiguration* GetMolecule() const override { return Molecule; }
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class MoleculeCounter : public G4VUserMoleculeCounter<MoleculeCounterIndex>
{
  public:
    MoleculeCounter(G4String = "MoleculeCounter");

    ~MoleculeCounter() override = default;

  public:
    void InitializeUser() override;

    std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*) const override;

    std::unique_ptr<G4VMoleculeCounterIndex> BuildIndex(const G4Track*,
                                                        const G4StepPoint*) const override;

    std::unique_ptr<G4VMoleculeCounterIndex>
    BuildSimpleIndex(const G4MolecularConfiguration*) const override;

    G4bool GetIgnoreMoleculePosition() const;

    void SetIgnoreMoleculePosition(G4bool);

    void SetNegativeCountsAreFatal(G4bool);

  protected:
    static G4ThreadLocal std::unique_ptr<G4Navigator> fNavigator;

    G4bool fIgnoreMoleculePosition{false};
    G4bool fIgnoreCopyNumbers{false};
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
