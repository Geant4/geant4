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
// $Id: G4RootNtupleManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for Root ntuples.
// It implements functions specific to Root ntuples.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4RootNtupleManager_h
#define G4RootNtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "globals.hh"

#include "tools/wroot/ntuple"

class G4RootMainNtupleManager;
class G4RootFileManager;

namespace tools {
namespace wroot {
class directory;
}
}

enum class G4NtupleCreateMode {
  kNoMergeBeforeOpen,
  kNoMergeAfterOpen,
  kMainBeforeOpen,
  kMainAfterOpen,
  kUndefined
};

// template specializations used by this class defined below

template <>
template <>
G4bool G4TNtupleManager<tools::wroot::ntuple>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value);


class G4RootNtupleManager : public G4TNtupleManager<tools::wroot::ntuple> 
{
  friend class G4RootAnalysisManager;
  friend class G4RootMainNtupleManager;

  public:
    explicit G4RootNtupleManager(const G4AnalysisManagerState& state, 
                                 G4int nofMainManagers = 0,
                                 G4bool rowWise = true);
    virtual ~G4RootNtupleManager();

   private:
    // Types alias
    using NtupleType = tools::wroot::ntuple;
    using NtupleDescriptionType = G4TNtupleDescription<NtupleType>;

    // Functions specific to the output type

    void SetNtupleDirectory(tools::wroot::directory* directory);
    void SetFileManager(std::shared_ptr<G4RootFileManager> fileManager);
    
    // Utility function  
    void CreateTNtuple(NtupleDescriptionType*  ntupleDescription);

    // Methods from the templated base class

    virtual void CreateTNtupleFromBooking(
                    NtupleDescriptionType*  ntupleDescription) final;

    virtual void FinishTNtuple(
                    NtupleDescriptionType*  ntupleDescription) final;

    virtual G4bool Reset(G4bool deleteNtuple);

    // Method for merging
    //
    G4bool Merge();

    // Access functions
    //
    const std::vector<NtupleDescriptionType*>& GetNtupleDescriptionVector() const;
    G4RootMainNtupleManager* GetMainNtupleManager(G4int index) const;
    unsigned int GetBasketSize() const;

    // Utility functions
    //
    void SetCreateMode();

    // data members
    //
    G4NtupleCreateMode fCreateMode;
    std::shared_ptr<G4RootFileManager> fFileManager;
    tools::wroot::directory*  fNtupleDirectory;
    std::vector<G4RootMainNtupleManager*>  fMainNtupleManagers;
};    

// inline functions

inline void 
G4RootNtupleManager::SetNtupleDirectory(tools::wroot::directory* directory) 
{ fNtupleDirectory = directory; }

inline void 
G4RootNtupleManager::SetFileManager(std::shared_ptr<G4RootFileManager> fileManager)
{ fFileManager = fileManager; }

inline const std::vector<G4TNtupleDescription<tools::wroot::ntuple>* >& 
G4RootNtupleManager::GetNtupleDescriptionVector() const
{ return fNtupleDescriptionVector; }

//_____________________________________________________________________________
template <>
template <>
inline G4bool G4TNtupleManager<tools::wroot::ntuple>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    //G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
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
  auto column = dynamic_cast<tools::wroot::ntuple::column_string* >(icolumn);
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

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "ntuple T column", description);
  }  
#endif
  return true;  
}

#endif


