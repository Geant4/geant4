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

// Manager class for CSV ntuples.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4CsvNtupleManager_h
#define G4CsvNtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "globals.hh"

#include "tools/wcsv_ntuple"

#include <memory>
#include <string_view>
#include <utility>

// Types alias
using CsvNtupleDescription = G4TNtupleDescription<tools::wcsv::ntuple, std::ofstream>;

class G4CsvFileManager;

class G4CsvNtupleManager : public G4TNtupleManager<tools::wcsv::ntuple,
                                                   std::ofstream>
{
  friend class G4CsvAnalysisManager;
  friend class G4CsvNtupleFileManager;

  public:
    explicit G4CsvNtupleManager(const G4AnalysisManagerState& state);
    G4CsvNtupleManager() = delete;
    ~G4CsvNtupleManager() override = default;

  private:
    // Functions specific to the output type

    // Set methods
    void SetFileManager(std::shared_ptr<G4CsvFileManager> fileManager);

    // Access to ntuple vector (needed for Write())
    const std::vector<CsvNtupleDescription*>& GetNtupleDescriptionVector() const;

    void SetIsCommentedHeader(G4bool isCommentedHeader);
    void SetIsHippoHeader(G4bool isHippoHeader);

    // Methods from the templated base class
    //
    void CreateTNtupleFromBooking(CsvNtupleDescription* ntupleDescription) final;

    void FinishTNtuple(CsvNtupleDescription* ntupleDescription, G4bool fromBooking) final;

    G4bool WriteHeader(tools::wcsv::ntuple* ntuple) const;

    // Static data members
    static constexpr std::string_view fkClass { "G4CsvNtupleManager" };

    // data members
    std::shared_ptr<G4CsvFileManager>  fFileManager { nullptr };
    G4bool  fIsCommentedHeader { true };
    G4bool  fIsHippoHeader { false };
};

// inline functions

inline void
G4CsvNtupleManager::SetFileManager(std::shared_ptr<G4CsvFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline const std::vector<G4TNtupleDescription<tools::wcsv::ntuple, std::ofstream>*>&
G4CsvNtupleManager::GetNtupleDescriptionVector() const
{ return fNtupleDescriptionVector; }

inline void G4CsvNtupleManager::SetIsCommentedHeader(G4bool isCommentedHeader)
{ fIsCommentedHeader = isCommentedHeader; }

inline void G4CsvNtupleManager::SetIsHippoHeader(G4bool isHippoHeader)
{ fIsHippoHeader = isHippoHeader; }

#endif
