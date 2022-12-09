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

// Manager class for Root ntuples.
// It implements functions specific to Root ntuples in sequential mode.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4RootNtupleManager_h
#define G4RootNtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "G4RootFileDef.hh"
#include "globals.hh"

#include "tools/wroot/ntuple"

#include <string_view>

class G4RootMainNtupleManager;
class G4NtupleBookingManager;
class G4RootFileManager;

// Types alias
using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;

// template specializations used by this class defined below

template <>
template <>
G4bool G4TNtupleManager<tools::wroot::ntuple, G4RootFile>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value);


class G4RootNtupleManager : public G4TNtupleManager<tools::wroot::ntuple,
                                                    G4RootFile>
{
  friend class G4RootAnalysisManager;
  friend class G4RootMainNtupleManager;
  friend class G4RootNtupleFileManager;
  friend class G4RootMpiNtupleFileManager;
  friend class G4RootMpiNtupleManager;

  public:
    G4RootNtupleManager(const G4AnalysisManagerState& state,
                const std::shared_ptr<G4NtupleBookingManager>& bookingManger,
                G4int nofMainManagers, G4int nofReducedFiles,
                G4bool rowWise, G4bool rowMode);
    G4RootNtupleManager() = delete;
    ~G4RootNtupleManager() override = default;

  private:
    // Methods from the templated base class
    void CreateTNtupleFromBooking(RootNtupleDescription* ntupleDescription) final;
    void FinishTNtuple(RootNtupleDescription* ntupleDescription, G4bool fromBooking) final;
    G4bool Reset() final;
    void Clear() final;
    virtual G4bool Merge();

    // Set functions
    void SetFileManager(const std::shared_ptr<G4RootFileManager>& fileManager);
    void SetNtupleFile(std::shared_ptr<G4RootFile> file);
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode);

    // New cycle option
    void SetNewCycle(G4bool value) final;

    // Access functions
    const std::vector<RootNtupleDescription*>& GetNtupleDescriptionVector() const;
    std::shared_ptr<G4RootMainNtupleManager> GetMainNtupleManager(G4int index) const;
    unsigned int GetBasketSize() const;
    unsigned int GetBasketEntries() const;

    // Static data members
    static constexpr std::string_view fkClass { "G4RootNtupleManager" };

    // Rata members
    std::shared_ptr<G4RootFileManager> fFileManager { nullptr };
    std::vector<std::shared_ptr<G4RootMainNtupleManager>>  fMainNtupleManagers;
    std::shared_ptr<G4RootFile> fNtupleFile { nullptr };
    G4bool fRowWise;
    G4bool fRowMode;
};

#include "G4RootNtupleManager.icc"

#endif


