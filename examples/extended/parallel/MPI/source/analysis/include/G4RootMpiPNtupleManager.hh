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

// Class for Root MPI pntuple management.
// The class handles ntuples on the processing MPI ranks when MPI ntuple merging 
// is activated.
// This class is temporarily provided with g4mpi,
// it will be integrated in Geant4 analysis category in future.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4RootMpiPNtupleManager_h
#define G4RootMpiPNtupleManager_h 1

#include "G4BaseNtupleManager.hh"
#include "G4RootMpiPNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include "tools/wroot/base_pntuple"

#include <vector>
#include <string_view>

using std::to_string;

class G4RootFileManager;

namespace tools {
namespace wroot {
class file;
class directory;
class impi_ntuple;
}
}

class G4RootMpiPNtupleManager : public G4BaseNtupleManager 
{
  friend class G4RootMpiNtupleFileManager;

  public:
    G4RootMpiPNtupleManager(const G4AnalysisManagerState& state, 
                            tools::impi* impi, G4int mpiRank, G4int destinationRank);
    ~G4RootMpiPNtupleManager();

  private:
    // Functions specific to the output type
    // void SetNtupleDirectory(tools::wroot::directory* directory);
    void SetFileManager(std::shared_ptr<G4RootFileManager> fileManager);

    // Methods to manipulate ntuples
    void CreateNtuple(G4RootMpiPNtupleDescription* ntupleDescription);
    void CreateNtuplesFromBooking(const std::vector<G4NtupleBooking*>& ntupleBookings);

    // Methods to create ntuples
    //
    virtual G4int CreateNtuple(G4NtupleBooking* booking) final;

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId (from base class)                    
    using G4BaseNtupleManager::FillNtupleIColumn;
    using G4BaseNtupleManager::FillNtupleFColumn;
    using G4BaseNtupleManager::FillNtupleDColumn;
    using G4BaseNtupleManager::FillNtupleSColumn;
    using G4BaseNtupleManager::AddNtupleRow;
    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)                      
    virtual G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) final;
    virtual G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) final;
    virtual G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) final;
    virtual G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId, 
                                     const G4String& value) final;
    virtual G4bool AddNtupleRow(G4int ntupleId) final;
    virtual G4bool Merge() final;

    // Clear all data
    virtual G4bool Reset() final;
    virtual void Clear() final;

    // Reset
    virtual G4bool Reset(G4bool deleteNtuple) final;

    // Activation option
    //
    virtual void  SetActivation(G4bool activation) final;
    virtual void  SetActivation(G4int ntupleId, G4bool activation) final;
    virtual G4bool  GetActivation(G4int ntupleId) const final;

    // New cycle option
    void SetNewCycle(G4bool value) final;
    G4bool GetNewCycle() const final;

    // Access methods
    virtual G4int GetNofNtuples() const final;

    // List ntuples
    G4bool List(std::ostream& output, G4bool onlyIfActive = true) final;

  private:
    G4RootMpiPNtupleDescription*  
      GetNtupleDescriptionInFunction(G4int id, G4String function, G4bool warn = true) const;
    tools::wroot::base_pntuple*  
      GetNtupleInFunction(G4int id, G4String function, G4bool warn = true) const;

    template <typename T> 
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);

    // Static data members
    static constexpr std::string_view fkClass { "G4RootMpiPNtupleManager" };

    // Data members
    // G4RootMpiMainNtupleManager* fMpiMainNtupleManager;
    std::shared_ptr<G4RootFileManager> fFileManager;
    std::vector<G4RootMpiPNtupleDescription*> fNtupleDescriptionVector;
    std::vector<tools::wroot::impi_ntuple*> fNtupleVector;
    tools::impi*  fImpi;
    G4int  fMpiRank;
    G4int  fDestinationRank;
    G4bool fNewCycle { false };
};

// inline functions

inline void 
G4RootMpiPNtupleManager::SetFileManager(std::shared_ptr<G4RootFileManager> fileManager)
{ fFileManager = fileManager; }

//_____________________________________________________________________________
template <>
inline G4bool G4RootMpiPNtupleManager::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
  } 

  if ( IsVerbose(G4Analysis::kVL4) ) {
    Message(G4Analysis::kVL4, "fill", "pntuple T column",
       " ntupleId " + to_string(ntupleId) + " columnId " + to_string(columnId) +
       " value " + G4Analysis::ToString(value));
  }

  auto ntuple = GetNtupleInFunction(ntupleId, "FillNtupleTColumn");
  if ( ! ntuple ) return false;

  auto index = columnId - fFirstNtupleColumnId;
  if ( index < 0 || index >= G4int(ntuple->columns().size()) ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::FillNtupleTColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  auto icolumn =  ntuple->columns()[index];
  auto column = dynamic_cast<tools::wroot::base_pntuple::column_string* >(icolumn);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << " Column type does not match: "
                << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    G4Exception("G4RootNtupleManager:FillNtupleColumn",
                "Analysis_W011", JustWarning, description);
    return false;
  } 

  column->fill(value);

  if ( IsVerbose(G4Analysis::kVL4) ) {
    Message(G4Analysis::kVL4, "done fill", "pntuple T column",
      " ntupleId " + to_string(ntupleId) +
      " columnId " + to_string(columnId) +
      " value " + value);
  }

  return true;  
}

//_____________________________________________________________________________
template <typename T>
G4bool G4RootMpiPNtupleManager::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const T& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
  }  

  if ( IsVerbose(G4Analysis::kVL4) ) {
    Message(G4Analysis::kVL4, "fill", "pntuple T column",
       " ntupleId " + to_string(ntupleId) +
       " columnId " + to_string(columnId) +
       " value " + G4Analysis::ToString(value));
  }

  // get ntuple
  auto ntuple = GetNtupleInFunction(ntupleId, "FillNtupleTColumn");
  if ( ! ntuple ) return false;

  // get generic column
  auto index = columnId - fFirstNtupleColumnId;
  if ( index < 0 || index >= G4int(ntuple->columns().size()) ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4TNtupleManager::FillNtupleTColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }
  auto icolumn =  ntuple->columns()[index];

  // get column and check its type
  auto column = dynamic_cast<tools::wroot::base_pntuple::column<T>* >(icolumn);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << " Column type does not match: "
                << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    G4Exception("G4TNtupleManager:FillNtupleTColumn",
                "Analysis_W011", JustWarning, description);
    return false;
  } 

  column->fill(value);

  if ( IsVerbose(G4Analysis::kVL4) ) {
    Message(G4Analysis::kVL4, "done fill", "pntuple T column",
       " ntupleId " + to_string(ntupleId) +
       " columnId " + to_string(columnId) +
       " value " + G4Analysis::ToString(value));
  }

  return true;  
}

#endif
