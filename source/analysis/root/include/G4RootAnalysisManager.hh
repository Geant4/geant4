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

#ifndef G4RootAnalysisManager_h
#define G4RootAnalysisManager_h 1

#include "G4ToolsAnalysisManager.hh"
#include "globals.hh"

#include "tools/wroot/ntuple"

#include <memory>
#include <string_view>

class G4RootAnalysisManager;
class G4RootFileManager;
class G4RootNtupleFileManager;
template <class T>
class G4ThreadLocalSingleton;

namespace tools {
namespace wroot {
class directory;
}
}

class G4RootAnalysisManager : public  G4ToolsAnalysisManager
{
  friend class G4RootMpiAnalysisManager;
  friend class G4ThreadLocalSingleton<G4RootAnalysisManager>;

  public:
    ~G4RootAnalysisManager() override;

    // Static methods
    static G4RootAnalysisManager* Instance();
    static G4bool IsInstance();

    // Access methods
    tools::wroot::ntuple* GetNtuple() const;
    tools::wroot::ntuple* GetNtuple(G4int ntupleId) const;

    // Iterators
    std::vector<tools::wroot::ntuple*>::iterator BeginNtuple();
    std::vector<tools::wroot::ntuple*>::iterator EndNtuple();
    std::vector<tools::wroot::ntuple*>::const_iterator BeginConstNtuple() const;
    std::vector<tools::wroot::ntuple*>::const_iterator EndConstNtuple() const;

    // MT/MPI
    void SetNtupleMerging(G4bool mergeNtuples, G4int nofReducedNtupleFiles = 0) override;
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode = true) override;
    void SetBasketSize(unsigned int basketSize) override;
    void SetBasketEntries(unsigned int basketEntries) override;

  private:
    G4RootAnalysisManager();

    // Static data members
    inline static G4ThreadLocal G4bool fgIsInstance { false };
    static constexpr std::string_view fkClass { "G4RootAnalysisManager" };

    // Data members
    std::shared_ptr<G4RootFileManager> fFileManager { nullptr };
    std::shared_ptr<G4RootNtupleFileManager> fNtupleFileManager { nullptr };
};

#include "G4RootAnalysisManager.icc"

#endif
