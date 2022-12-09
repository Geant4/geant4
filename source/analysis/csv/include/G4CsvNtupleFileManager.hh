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

// Manager class for ntuple Csv file output.
//
// Author: Ivana Hrivnacova, 15/09/2020 (ivana@ipno.in2p3.fr)

#ifndef G4CsvNtupleFileManager_h
#define G4CsvNtupleFileManager_h 1

#include "G4VNtupleFileManager.hh"
#include "globals.hh"

#include <string_view>
#include <utility>

class G4CsvFileManager;
class G4CsvNtupleManager;
class G4VNtupleManager;
class G4NtupleBookingManager;

class G4CsvNtupleFileManager : public G4VNtupleFileManager
{
  friend class G4CsvAnalysisManager;

  public:
    explicit G4CsvNtupleFileManager(const G4AnalysisManagerState& state);
    G4CsvNtupleFileManager() = delete;
    ~G4CsvNtupleFileManager() override = default;

    std::shared_ptr<G4VNtupleManager> CreateNtupleManager() override;

    // Methods to be performed at file management
    G4bool ActionAtOpenFile(const G4String& fileName) override;
    G4bool ActionAtWrite() override;
    G4bool ActionAtCloseFile() override;
    G4bool Reset() override;

    void SetFileManager(std::shared_ptr<G4CsvFileManager> fileManager);

    std::shared_ptr<G4CsvNtupleManager> GetNtupleManager() const;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4CsvNtupleFileManager" };

    // Data members
    std::shared_ptr<G4CsvFileManager> fFileManager { nullptr };
    std::shared_ptr<G4CsvNtupleManager> fNtupleManager { nullptr };
};

// inline functions

inline void G4CsvNtupleFileManager::SetFileManager(
  std::shared_ptr<G4CsvFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline std::shared_ptr<G4CsvNtupleManager> G4CsvNtupleFileManager::GetNtupleManager() const
{ return fNtupleManager; }

#endif

