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
// $Id:$

// Class for Root pntuple management.
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#ifndef G4RootPNtupleManager_h
#define G4RootPNtupleManager_h 1

#include "G4BaseNtupleManager.hh"
#include "G4RootPNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AutoLock.hh"
#include "globals.hh"

#include <vector>

class G4RootMainNtupleManager;

namespace tools {
namespace wroot {
class file;
class ntuple;
class imt_ntuple;
}
}

// Mutex implementation as in pwroot.cpp
// with replacement tools::mutex -> G4Mutex

class mutex : public virtual tools::wroot::imutex {
  typedef tools::wroot::imutex parent;
public:
  virtual bool lock() {
    // G4cout << "!!! Mutex lock" << G4endl;
    m_mutex.lock(); 
    return true;}
  virtual bool unlock() {
    m_mutex.unlock(); 
    // G4cout << "!!! Mutex unlock" << G4endl;
    return true; }
  //virtual bool trylock() {return m_mutex.trylock();}
public:
  mutex(G4AutoLock& a_mutex):m_mutex(a_mutex){}
  virtual ~mutex(){}
protected:
  mutex(const mutex& a_from):parent(a_from),m_mutex(a_from.m_mutex){}
  mutex& operator=(const mutex&){return *this;}
protected:
  G4AutoLock& m_mutex;
};

class G4RootPNtupleManager : public G4BaseNtupleManager 
{
  friend class G4RootAnalysisManager;

  public:
    explicit G4RootPNtupleManager(G4RootMainNtupleManager* main,
                                  const G4AnalysisManagerState& state);
    ~G4RootPNtupleManager();

  private:
    enum class G4PNtupleCreateMode {
      kSlaveBeforeOpen,
      kSlaveAfterOpen,
      kUndefined
    };

    // Methods to manipulate ntuples
    void CreateNtuple(G4RootPNtupleDescription* ntupleDescription,
                      tools::wroot::ntuple* mainNtuple);
    void CreateNtuplesFromMain();

    // Methods to create ntuples
    //
    virtual G4int CreateNtuple(const G4String& name, const G4String& title) final;

    // Create columns in the last created ntuple (from base class)
    using G4BaseNtupleManager::CreateNtupleIColumn;
    using G4BaseNtupleManager::CreateNtupleFColumn;
    using G4BaseNtupleManager::CreateNtupleDColumn;
    using G4BaseNtupleManager::CreateNtupleSColumn;
    using G4BaseNtupleManager::FinishNtuple; 

    // Create columns in the ntuple with given id
    virtual G4int CreateNtupleIColumn(G4int ntupleId, 
                    const G4String& name, std::vector<int>* vector) override;
    virtual G4int CreateNtupleFColumn(G4int ntupleId, 
                    const G4String& name, std::vector<float>* vector) override;
    virtual G4int CreateNtupleDColumn(G4int ntupleId, 
                    const G4String& name, std::vector<double>* vector) override;
    virtual G4int CreateNtupleSColumn(G4int ntupleId, const G4String& name) override;
    virtual void  FinishNtuple(G4int ntupleId) override;   

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

    // Reset
    virtual G4bool Reset(G4bool deleteNtuple) final;

    // Activation option
    //
    virtual void  SetActivation(G4bool activation) final;
    virtual void  SetActivation(G4int ntupleId, G4bool activation) final;
    virtual G4bool  GetActivation(G4int ntupleId) const final;
    virtual G4bool  IsEmpty() const final;

    // Access methods
    virtual G4int GetNofNtuples() const final;
    virtual G4int GetNofNtupleBookings() const final;

  private:
    G4RootPNtupleDescription*  
      GetNtupleDescriptionInFunction(G4int id, G4String function, G4bool warn = true) const;
    tools::wroot::base_pntuple*  
      GetNtupleInFunction(G4int id, G4String function, G4bool warn = true) const;
    tools::wroot::ntuple*  
      GetMainNtupleInFunction(G4int id, G4String function, G4bool warn = true) const;

    template <typename T> 
    G4int CreateNtupleTColumn(G4int ntupleId, 
                    const G4String& name, std::vector<T>* vector);

    template <typename T> 
    G4int CreateNtupleTColumn(
                    const G4String& name, std::vector<T>* vector);

    template <typename T> 
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);

    template <typename T> 
    G4bool FillNtupleTColumn(G4int columnId, const T& value);

    // Data members
    G4PNtupleCreateMode  fCreateMode;
    G4RootMainNtupleManager*  fMainNtupleManager;
    std::vector<G4RootPNtupleDescription*> fNtupleDescriptionVector;
    std::vector<tools::wroot::imt_ntuple*> fNtupleVector;
};

// inline functions

//_____________________________________________________________________________
template <typename T>
G4int G4RootPNtupleManager::CreateNtupleTColumn(
  G4int ntupleId, const G4String& name, std::vector<T>* vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "pntuple T column", description);
  }  
#endif
  
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "CreateNtupleTColumn");
  if ( ! ntupleDescription )  return G4Analysis::kInvalidId;   
    
  // Save column info in booking
  auto& ntupleBooking = ntupleDescription->fNtupleBooking;
  auto index = ntupleBooking.columns().size();
  if ( ! vector )
    ntupleBooking.template add_column<T>(name);
  else
    ntupleBooking.template add_column<T>(name, *vector);
 
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL2()->Message("create", "pntuple T column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       

}

//_____________________________________________________________________________
template <typename T> 
G4int G4RootPNtupleManager::CreateNtupleTColumn(
  const G4String& name, std::vector<T>* vector)
{
  auto ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleTColumn<T>(ntupleId, name, vector);
}

//_____________________________________________________________________________
template <>
inline G4bool G4RootPNtupleManager::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const std::string& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
  } 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "pntuple T column", description);
  }  
#endif

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

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("done fill", "pntuple T column", description);
  }  
#endif
  return true;  
}

//_____________________________________________________________________________
template <typename T>
G4bool G4RootPNtupleManager::FillNtupleTColumn(
  G4int ntupleId, G4int columnId, const T& value)
{
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    G4cout << "Skipping FillNtupleIColumn for " << ntupleId << G4endl; 
    return false; 
  }  

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "pntuple T column", description);
  }  
#endif

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

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("done fill", "pntuple T column", description);
  }  
#endif
  return true;  
}

//_____________________________________________________________________________
template <typename T> 
G4bool G4RootPNtupleManager::FillNtupleTColumn(G4int columnId, const T& value) {
  return FillNtupleTColumn(0, columnId, value);
}

#endif

