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
  friend class G4RootNtupleFileManager;

  public:
    explicit G4RootPNtupleManager(const G4AnalysisManagerState& state,
                                  std::shared_ptr<G4NtupleBookingManager> bookingManger,
                                  std::shared_ptr<G4RootMainNtupleManager> main,
                                  G4bool rowWise, G4bool rowMode);
    ~G4RootPNtupleManager();

  private:
    // Methods to manipulate ntuples
    void CreateNtupleFromMain(G4RootPNtupleDescription* ntupleDescription,
                              tools::wroot::ntuple* mainNtuple);
    // void CreateNtupleFromMain(G4NtupleBooking* g4NtupleBooking,
    //                           tools::wroot::ntuple* mainNtuple);
    void CreateNtuplesFromMain();

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

    // Reset
    virtual G4bool Reset(G4bool deleteNtuple) final;

    // Activation option
    //
    virtual void  SetActivation(G4bool activation) final;
    virtual void  SetActivation(G4int ntupleId, G4bool activation) final;
    virtual G4bool  GetActivation(G4int ntupleId) const final;

    // Access methods
    virtual G4int GetNofNtuples() const final;

    // Set methods
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode);

  private:
    G4RootPNtupleDescription*  
      GetNtupleDescriptionInFunction(G4int id, G4String function, G4bool warn = true) const;
    tools::wroot::base_pntuple*  
      GetNtupleInFunction(G4int id, G4String function, G4bool warn = true) const;
    tools::wroot::ntuple*  
      GetMainNtupleInFunction(G4int id, G4String function, G4bool warn = true) const;

    template <typename T> 
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);

    template <typename T> 
    G4bool FillNtupleTColumn(G4int columnId, const T& value);

    // Data members
    std::shared_ptr<G4NtupleBookingManager> fBookingManager;
    std::shared_ptr<G4RootMainNtupleManager>  fMainNtupleManager;
    std::vector<G4RootPNtupleDescription*> fNtupleDescriptionVector;
    std::vector<tools::wroot::imt_ntuple*> fNtupleVector;
    G4bool fRowWise;
    G4bool fRowMode;
    G4bool fCreateNtuples;
};

#include "G4RootPNtupleManager.icc"

#endif
