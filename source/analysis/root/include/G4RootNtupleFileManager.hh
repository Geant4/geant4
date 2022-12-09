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

// Manager class for ntuple Root file output,
// provides handling Root ntuple managers in parallel processing.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4RootNtupleFileManager_h
#define G4RootNtupleFileManager_h 1

#include "G4VNtupleFileManager.hh"

#include <string_view>
#include <utility>

class G4RootFileManager;
class G4RootNtupleManager;
class G4RootPNtupleManager;
class G4VNtupleManager;
class G4NtupleBookingManager;

class G4RootNtupleFileManager : public G4VNtupleFileManager
{
  friend class G4RootMpiNtupleFileManager;

  public:
    explicit G4RootNtupleFileManager(const G4AnalysisManagerState& state);
    G4RootNtupleFileManager() = delete;
    ~G4RootNtupleFileManager() override;

    // MT/MPI
    void SetNtupleMerging(G4bool mergeNtuples, G4int nofReducedNtupleFiles = 0) override;
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true) override;
    void SetBasketSize(unsigned int basketSize) override;
    void SetBasketEntries(unsigned int basketEntries) override;

    // virtual methods from base class
    G4bool ActionAtOpenFile(const G4String& fileName) override;
    G4bool ActionAtWrite() override;
    G4bool ActionAtCloseFile() override;
    G4bool Reset() override;
    G4bool IsNtupleMergingSupported() const override;

    std::shared_ptr<G4VNtupleManager> CreateNtupleManager() override;

    void SetFileManager(std::shared_ptr<G4RootFileManager> fileManager);

    G4NtupleMergeMode GetMergeMode() const override;
    std::shared_ptr<G4RootNtupleManager> GetNtupleManager() const;

  private:
    // Static data members
    static G4RootNtupleFileManager* fgMasterInstance;

    void  SetNtupleMergingMode(G4bool mergeNtuples, G4int nofNtupleFiles);
    G4int GetNtupleFileNumber();
    G4bool CloseNtupleFiles();

    // Static data members
    static constexpr std::string_view fkClass { "G4RootNtupleFileManager" };

    // data members
    G4bool  fIsInitialized { false };
    G4int   fNofNtupleFiles { 0 };
    G4bool  fNtupleRowWise { false };
    G4bool  fNtupleRowMode { true };
    G4NtupleMergeMode  fNtupleMergeMode { G4NtupleMergeMode::kNone };
    std::shared_ptr<G4RootNtupleManager>  fNtupleManager { nullptr };
    std::shared_ptr<G4RootPNtupleManager> fSlaveNtupleManager { nullptr };
    std::shared_ptr<G4RootFileManager>    fFileManager { nullptr };
};

inline void G4RootNtupleFileManager::SetFileManager(
  std::shared_ptr<G4RootFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline G4NtupleMergeMode G4RootNtupleFileManager::GetMergeMode() const
{ return fNtupleMergeMode; }

inline G4bool G4RootNtupleFileManager::IsNtupleMergingSupported() const
{ return true; }

inline std::shared_ptr<G4RootNtupleManager> G4RootNtupleFileManager::GetNtupleManager() const
{ return fNtupleManager; }

#endif
