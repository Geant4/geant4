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
// $Id$

// Manager class for Hdf5 ntuples 
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5NtupleManager_h
#define G4Hdf5NtupleManager_h 1

#include "G4TNtupleManager.hh"
#include "globals.hh"

#include "tools/hdf5/ntuple"

#include <vector>
#include <memory>

class G4Hdf5FileManager;

// template specialization used by this class defined below

template <>
template <>
G4bool G4TNtupleManager<tools::hdf5::ntuple>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value);


class G4Hdf5NtupleManager : public G4TNtupleManager<tools::hdf5::ntuple> 
{
  friend class G4Hdf5AnalysisManager;

  public:
    explicit G4Hdf5NtupleManager(const G4AnalysisManagerState& state);
    ~G4Hdf5NtupleManager();

  private:
    // Types alias
    using NtupleType = tools::hdf5::ntuple;
    using NtupleDescriptionType = G4TNtupleDescription<NtupleType>;

    // Set methods
    void SetFileManager(std::shared_ptr<G4Hdf5FileManager> fileManager);
    
    // Access to ntuple vector (needed for Write())
    const std::vector<NtupleDescriptionType*>& GetNtupleDescriptionVector() const;

    // Methods from the templated base class
    //
    virtual void CreateTNtuple(
                    NtupleDescriptionType*  ntupleDescription,
                    const G4String& name, const G4String& title) final;
    virtual void CreateTNtupleFromBooking(
                    NtupleDescriptionType* ntupleDescription) final;

    virtual void FinishTNtuple(
                    NtupleDescriptionType* ntupleDescription) final;

    // data members
    //
    std::shared_ptr<G4Hdf5FileManager>  fFileManager;
};

// inline functions

inline void 
G4Hdf5NtupleManager::SetFileManager(std::shared_ptr<G4Hdf5FileManager> fileManager)
{ fFileManager = fileManager; }

inline const std::vector<G4TNtupleDescription<tools::hdf5::ntuple>*>& 
G4Hdf5NtupleManager::GetNtupleDescriptionVector() const
{ return fNtupleDescriptionVector; }

template <>
template <>
inline G4bool G4TNtupleManager<tools::hdf5::ntuple>::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    //G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
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
  auto column = dynamic_cast<tools::hdf5::ntuple::column_string* >(icolumn);
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

