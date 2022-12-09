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

// Class for Root main ntuple management.
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#ifndef G4RootMainNtupleManager_h
#define G4RootMainNtupleManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "G4NtupleBookingManager.hh"
#include "G4RootFileDef.hh"
#include "G4RootNtupleManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include <string_view>
#include <utility>
#include <vector>

class G4RootFileManager;
class G4RootNtupleManager;
class G4NtupleBookingManager;

namespace tools {
namespace wroot {
class ntuple;
}
}

// Types alias
using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;
using RootMainNtupleDescription = std::pair<RootNtupleDescription*, std::shared_ptr<G4RootFile>>;

class G4RootMainNtupleManager : public G4BaseAnalysisManager
{
  friend class G4RootPNtupleManager;
  friend class G4RootNtupleManager;

  public:
    G4RootMainNtupleManager(
               G4RootNtupleManager* ntupleBuilder,
               std::shared_ptr<G4NtupleBookingManager> bookingManager,
               G4bool rowWise, G4int fileNumber,
               const G4AnalysisManagerState& state);
    G4RootMainNtupleManager() = delete;
    ~G4RootMainNtupleManager() override = default;

    void CreateNtuplesFromBooking();

  protected:
    // Methods to manipulate ntuples
    void   CreateNtuple(RootNtupleDescription* ntupleDescription, G4bool warn = true);
    G4bool Merge();
    G4bool Reset();
    void ClearData();

    // Set/get methods
    void SetFileManager(std::shared_ptr<G4RootFileManager> fileManager);
    void SetRowWise(G4bool rowWise);
    std::shared_ptr<G4RootFile> GetNtupleFile(RootNtupleDescription* ntupleDescription) const;

    // New cycle option
    void SetNewCycle(G4bool value);
    G4bool GetNewCycle() const;

    // Access functions
    const std::vector<tools::wroot::ntuple*>& GetNtupleVector()
      { return fNtupleVector; }
    unsigned int GetBasketEntries() const;

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4RootMainNtupleManager" };

    // Data members
    G4RootNtupleManager* fNtupleBuilder { nullptr };
    std::shared_ptr<G4NtupleBookingManager> fBookingManager { nullptr };
    std::shared_ptr<G4RootFileManager> fFileManager { nullptr };
    G4bool  fRowWise;
    G4int  fFileNumber;
    std::vector<tools::wroot::ntuple*>  fNtupleVector;
    std::vector<RootMainNtupleDescription> fNtupleDescriptionVector;
    G4bool fNewCycle { false };
};

inline void  G4RootMainNtupleManager::SetFileManager(
  std::shared_ptr<G4RootFileManager> fileManager)
{
  fFileManager = std::move(fileManager);
}

inline void G4RootMainNtupleManager::SetRowWise(G4bool rowWise)
{ fRowWise = rowWise; }

inline unsigned int G4RootMainNtupleManager::GetBasketEntries() const
{ return fNtupleBuilder->GetBasketEntries(); }

inline void G4RootMainNtupleManager::SetNewCycle(G4bool value)
{ fNewCycle = value; }

inline G4bool G4RootMainNtupleManager::GetNewCycle() const
{ return fNewCycle; }

#endif
