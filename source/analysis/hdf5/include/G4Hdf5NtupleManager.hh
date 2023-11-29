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

// Manager class for Hdf5 ntuples
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5NtupleManager_h
#define G4Hdf5NtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include "toolx/hdf5/ntuple"

#include <memory>
#include <utility>
#include <vector>

class G4Hdf5FileManager;

// Types alias
using G4Hdf5File = std::tuple<hid_t, hid_t, hid_t>;
using Hdf5NtupleDescription = G4TNtupleDescription<toolx::hdf5::ntuple, G4Hdf5File>;

// template specialization used by this class defined below

template <>
template <>
G4bool G4TNtupleManager<toolx::hdf5::ntuple, G4Hdf5File>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value);


class G4Hdf5NtupleManager : public G4TNtupleManager<toolx::hdf5::ntuple,
                                                    G4Hdf5File>
{
  friend class G4Hdf5AnalysisManager;
  friend class G4Hdf5NtupleFileManager;

  public:
    explicit G4Hdf5NtupleManager(const G4AnalysisManagerState& state);
    G4Hdf5NtupleManager() = delete;
    ~G4Hdf5NtupleManager() override = default;

  private:
    // Set methods
    void SetFileManager(std::shared_ptr<G4Hdf5FileManager> fileManager);

    // Access to ntuple vector (needed for Write())
    const std::vector<Hdf5NtupleDescription*>& GetNtupleDescriptionVector() const;

    // Utility function
    void CreateTNtuple(Hdf5NtupleDescription* ntupleDescription, G4bool warn);

    // Methods from the templated base class
    //
    void CreateTNtupleFromBooking(Hdf5NtupleDescription* ntupleDescription) final;

    void FinishTNtuple(Hdf5NtupleDescription* ntupleDescription, G4bool fromBooking) final;

    // Static data members
    static constexpr std::string_view fkClass { "G4Hdf5NtupleManager" };

    // Data members
    std::shared_ptr<G4Hdf5FileManager>  fFileManager { nullptr };
};

// inline functions

using std::to_string;

inline void
G4Hdf5NtupleManager::SetFileManager(std::shared_ptr<G4Hdf5FileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline const std::vector<Hdf5NtupleDescription*>&
G4Hdf5NtupleManager::GetNtupleDescriptionVector() const
{ return fNtupleDescriptionVector; }

template <>
template <>
inline G4bool G4TNtupleManager<toolx::hdf5::ntuple, G4Hdf5File>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    //G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl;
    return false;
  }

  // get ntuple
  auto ntuple = GetNtupleInFunction(ntupleId, "FillNtupleTColumn");
  if (ntuple == nullptr) return false;

  // get generic column
  auto index = columnId - fFirstNtupleColumnId;
  if ( index < 0 || index >= G4int(ntuple->columns().size()) ) {
    G4Analysis::Warn(
      "ntupleId " + to_string(ntupleId) + " columnId " + to_string(columnId) +
      " does not exist.", fkClass, "FillNtupleTColumn");
    return false;
  }
  auto icolumn =  ntuple->columns()[index];

  // get column and check its type
  auto column = dynamic_cast<toolx::hdf5::ntuple::column_string* >(icolumn);
  if (column == nullptr) {
    G4Analysis::Warn(
      "Column type does not match: ntupleId " + to_string(ntupleId) +
      " columnId " + to_string(columnId) + " value " + value,
      fkClass, "FillNtupleTColumn");
    return false;
  }

  column->fill(value);

  if ( IsVerbose(G4Analysis::kVL4) ) {
    Message(G4Analysis::kVL4, "fill", "ntuple T column",
      " ntupleId " + to_string(ntupleId) +
      " columnId " + to_string(columnId) +
      " value " + G4Analysis::ToString(value));
  }

  return true;
}

#endif
