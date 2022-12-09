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
#include "G4RootMainNtupleManager.hh"
#include "G4RootPNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AutoLock.hh"
#include "globals.hh"

#include <vector>
#include <string_view>

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
  using parent = tools::wroot::imutex;
public:
  bool lock() override
  {
    // G4cout << "!!! Mutex lock" << G4endl;
    m_mutex.lock();
    return true;
  }
  bool unlock() override
  {
    m_mutex.unlock();
    // G4cout << "!!! Mutex unlock" << G4endl;
    return true;
  }
  //virtual bool trylock() {return m_mutex.trylock();}
public:
  mutex(G4AutoLock& a_mutex):m_mutex(a_mutex){}
  ~mutex() override = default;

protected:
  mutex(const mutex& a_from) = default;
  mutex& operator=(const mutex&){return *this;}
protected:
  G4AutoLock& m_mutex;
};

class G4RootPNtupleManager : public G4BaseNtupleManager
{
  friend class G4RootAnalysisManager;
  friend class G4RootNtupleFileManager;

  public:
    G4RootPNtupleManager(const G4AnalysisManagerState& state,
                         std::shared_ptr<G4NtupleBookingManager> bookingManger,
                         std::shared_ptr<G4RootMainNtupleManager> main,
                         G4bool rowWise, G4bool rowMode);
    G4RootPNtupleManager() = delete;
    ~G4RootPNtupleManager() override;

  private:
    // Methods to manipulate ntuples
    void CreateNtupleFromMain(G4RootPNtupleDescription* ntupleDescription,
                              tools::wroot::ntuple* mainNtuple);
    void CreateNtupleDescriptionsFromBooking();
    void CreateNtuplesFromMain();
    void CreateNtuplesIfNeeded();

    // Methods to create ntuples
    //
    G4int CreateNtuple(G4NtupleBooking* booking) final;

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId (from base class)
    using G4BaseNtupleManager::FillNtupleIColumn;
    using G4BaseNtupleManager::FillNtupleFColumn;
    using G4BaseNtupleManager::FillNtupleDColumn;
    using G4BaseNtupleManager::FillNtupleSColumn;
    using G4BaseNtupleManager::AddNtupleRow;
    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)
    G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) final;
    G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) final;
    G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) final;
    G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId, const G4String& value) final;
    G4bool AddNtupleRow(G4int ntupleId) final;
    virtual G4bool Merge() final;

    virtual G4bool Reset();
    void Clear() final;

    // Activation option
    //
    void SetActivation(G4bool activation) final;
    void SetActivation(G4int ntupleId, G4bool activation) final;
    G4bool GetActivation(G4int ntupleId) const final;

    // New cycle option
    void SetNewCycle(G4bool value) final;
    G4bool GetNewCycle() const final;

    // Access methods
    G4int GetNofNtuples() const final;

    // Set methods
    void SetNtupleRowWise(G4bool rowWise, G4bool rowMode);

    // List ntuples
    G4bool List(std::ostream& output, G4bool onlyIfActive = true) final;

  private:
    G4RootPNtupleDescription*
      GetNtupleDescriptionInFunction(G4int id, std::string_view function, G4bool warn = true) const;
    tools::wroot::base_pntuple*
      GetNtupleInFunction(G4int id, std::string_view function, G4bool warn = true) const;
    tools::wroot::ntuple*
      GetMainNtupleInFunction(G4int id, std::string_view function, G4bool warn = true) const;

    template <typename T>
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);

    template <typename T>
    G4bool FillNtupleTColumn(G4int columnId, const T& value);

    // Static data members
    static constexpr std::string_view fkClass { "G4RootPNtupleManager" };

    // Data members
    std::shared_ptr<G4NtupleBookingManager> fBookingManager;
    std::shared_ptr<G4RootMainNtupleManager>  fMainNtupleManager;
    std::vector<G4RootPNtupleDescription*> fNtupleDescriptionVector;
    std::vector<tools::wroot::imt_ntuple*> fNtupleVector;
    G4bool fRowWise;
    G4bool fRowMode;
    G4bool fCreateNtuples { true };
    G4bool fNewCycle { false };
};

#include "G4RootPNtupleManager.icc"

#endif
