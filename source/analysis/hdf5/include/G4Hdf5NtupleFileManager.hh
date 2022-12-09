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

// Manager class for ntuple Hdf5 file output.
//
// Author: Ivana Hrivnacova, 15/09/2020 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5NtupleFileManager_h
#define G4Hdf5NtupleFileManager_h 1

#include "G4VNtupleFileManager.hh"
#include "globals.hh"

#include <string_view>
#include <utility>

class G4Hdf5FileManager;
class G4Hdf5NtupleManager;
class G4VNtupleManager;
class G4NtupleBookingManager;

class G4Hdf5NtupleFileManager : public G4VNtupleFileManager
{
  friend class G4Hdf5AnalysisManager;

  public:
    explicit G4Hdf5NtupleFileManager(const G4AnalysisManagerState& state);
    G4Hdf5NtupleFileManager() = delete;
    ~G4Hdf5NtupleFileManager() override = default;

    std::shared_ptr<G4VNtupleManager> CreateNtupleManager() override;

    // Methods to be performed at file management
    G4bool ActionAtOpenFile(const G4String& fileName) override;
    G4bool ActionAtWrite() override;
    G4bool ActionAtCloseFile() override;
    G4bool Reset() override;

    void SetFileManager(std::shared_ptr<G4Hdf5FileManager> fileManager);

    std::shared_ptr<G4Hdf5NtupleManager> GetNtupleManager() const;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4Hdf5NtupleFileManager" };

    // Data members
    std::shared_ptr<G4Hdf5FileManager> fFileManager { nullptr };
    std::shared_ptr<G4Hdf5NtupleManager> fNtupleManager { nullptr };
};

inline void G4Hdf5NtupleFileManager::SetFileManager(
  std::shared_ptr<G4Hdf5FileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline std::shared_ptr<G4Hdf5NtupleManager> G4Hdf5NtupleFileManager::GetNtupleManager() const
{ return fNtupleManager; }

#endif
