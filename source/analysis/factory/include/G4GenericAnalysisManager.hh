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

// The main manager for Root analysis.
// It delegates most of functions to the object specific managers.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4GenericAnalysisManager_h
#define G4GenericAnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "G4THnManager.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include <memory>
#include <string_view>

class G4GenericAnalysisManager;
class G4GenericAnalysisMessenger;
class G4GenericFileManager;
class G4HnInformation;
class G4VNtupleFileManager;
template <class T>
class G4ThreadLocalSingleton;

class G4GenericAnalysisManager : public  G4ToolsAnalysisManager
{
  friend class G4RootMpiAnalysisManager;
  friend class G4ThreadLocalSingleton<G4GenericAnalysisManager>;

  public:
    ~G4GenericAnalysisManager() override;

    // Static methods
    static G4GenericAnalysisManager* Instance();
    static G4bool IsInstance();

    // MT/MPI
    void SetNtupleMerging(G4bool mergeNtuples, G4int nofReducedNtupleFiles = 0) override;
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true) override;
    void SetBasketSize(unsigned int basketSize) override;
    void SetBasketEntries(unsigned int basketEntries) override;

    // write in an extra file
    G4bool WriteH1(G4int id, const G4String& fileName);
    G4bool WriteH2(G4int id, const G4String& fileName);
    G4bool WriteH3(G4int id, const G4String& fileName);
    G4bool WriteP1(G4int id, const G4String& fileName);
    G4bool WriteP2(G4int id, const G4String& fileName);

    // Set default output type (backward compatibility)
    // this type will be used for file names without extension
    void SetDefaultFileType(const G4String& value);
    G4String GetDefaultFileType() const;

  protected:
    // Virtual methods from base class
    G4bool OpenFileImpl(const G4String& fileName) override;
    // File manager access
    std::shared_ptr<G4VFileManager> GetFileManager(const G4String& fileName) final;

  private:
    G4GenericAnalysisManager();
    // Methods
    void CreateNtupleFileManager(const G4String& fileName);

    // Static data members
    inline static G4GenericAnalysisManager* fgMasterInstance { nullptr };
    inline static G4ThreadLocal G4bool fgIsInstance { false };
    static constexpr std::string_view fkClass { "G4GenericAnalysisManager" };

    // Data members
    std::unique_ptr<G4GenericAnalysisMessenger>  fMessenger;
    std::shared_ptr<G4GenericFileManager> fFileManager { nullptr };
    // add G4GenericNtupleManager
    // this class will be analogic to file managers but with ntuples
    std::shared_ptr<G4VNtupleFileManager> fNtupleFileManager { nullptr };

    // Data members
    G4bool  fIsNtupleMergingSet { false };
    G4int   fNofNtupleFiles { 0 };
    G4bool  fMergeNtuples { false };
    G4bool  fNtupleRowWise { false };
    G4bool  fNtupleRowMode { true };
    unsigned int fBasketSize { G4Analysis::kDefaultBasketSize };
    unsigned int fBasketEntries { G4Analysis::kDefaultBasketEntries };
};

#include "G4GenericAnalysisManager.icc"

#endif
