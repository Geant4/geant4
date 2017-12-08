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
// $Id: G4CsvNtupleManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Class template for ntuple managers for all output types.
//
// Author: Ivana Hrivnacova, 19/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4TNtupleManager_h
#define G4TNtupleManager_h 1

#include "G4BaseNtupleManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include <vector>

template <typename TNTUPLE>
class G4TNtupleManager : public G4BaseNtupleManager {

  public:
    explicit G4TNtupleManager(const G4AnalysisManagerState& state);
    ~G4TNtupleManager();

  protected:
    // Methods to manipulate ntuples  
    void CreateNtuplesFromBooking();  
    G4bool IsEmpty() const;
    G4bool Reset(G4bool deleteNtuple);

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
                    const G4String& name, std::vector<int>* vector) final;
    virtual G4int CreateNtupleFColumn(G4int ntupleId, 
                    const G4String& name, std::vector<float>* vector) final;
    virtual G4int CreateNtupleDColumn(G4int ntupleId, 
                    const G4String& name, std::vector<double>* vector) final;
    virtual G4int CreateNtupleSColumn(G4int ntupleId, const G4String& name) final;
    virtual void  FinishNtuple(G4int ntupleId) final;   

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

    // Activation option
    //
    virtual void  SetActivation(G4bool activation) final;
    virtual void  SetActivation(G4int ntupleId, G4bool activation) final;
    virtual G4bool  GetActivation(G4int ntupleId) const final;

    // Access methods
    TNTUPLE* GetNtuple() const;
    TNTUPLE* GetNtuple(G4int ntupleId) const;
    virtual G4int GetNofNtuples() const final;
    virtual G4int GetNofNtupleBookings() const override;

    // Iterators
    typename std::vector<TNTUPLE*>::iterator BeginNtuple();  
    typename std::vector<TNTUPLE*>::iterator EndNtuple();
    typename std::vector<TNTUPLE*>::const_iterator BeginConstNtuple() const;
    typename std::vector<TNTUPLE*>::const_iterator EndConstNtuple() const;
 
    // Data members
    std::vector<G4TNtupleDescription<TNTUPLE>*> fNtupleDescriptionVector;
    std::vector<TNTUPLE*> fNtupleVector;

  private:
    // methods
   
    // Fuctions which are specific to output type
    //
    virtual void CreateTNtupleFromBooking(
                    G4TNtupleDescription<TNTUPLE>* ntupleDescription) = 0;

    virtual void FinishTNtuple(
                    G4TNtupleDescription<TNTUPLE>* ntupleDescription) = 0;
    
    void FinishTNtupleNew(G4TNtupleDescription<TNTUPLE>* ntupleDescription);

    // Common implementation
    //

    G4TNtupleDescription<TNTUPLE>*  GetNtupleDescriptionInFunction(G4int id, 
                                        G4String function,
                                        G4bool warn = true) const;
    TNTUPLE*  GetNtupleInFunction(G4int id, 
                                  G4String function,
                                  G4bool warn = true) const;

    // template functions for creating/filling ntuple columns

    template <typename T> 
    G4int CreateNtupleTColumn(G4int ntupleId, 
                    const G4String& name, std::vector<T>* vector);

    template <typename T> 
    G4bool FillNtupleTColumn(G4int ntupleId, G4int columnId, const T& value);
};

#include "G4TNtupleManager.icc"

#endif

