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

// Class template for read ntuple managers for all output types.
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4TRNtupleManager_h
#define G4TRNtupleManager_h 1

#include "G4BaseRNtupleManager.hh"
#include "G4TRNtupleDescription.hh"
#include "globals.hh"

#include <vector>

template <typename TNTUPLE>
class G4TRNtupleManager : public G4BaseRNtupleManager
{
  protected:
    explicit G4TRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4TRNtupleManager();

    // Methods to manipulate ntuples  
    G4bool IsEmpty() const;
    G4bool Reset();

    // Access methods
    TNTUPLE* GetNtuple() const;
    TNTUPLE* GetNtuple(G4int ntupleId) const;

    // Functions independent from the output type 
    //
    // Methods to read ntuple from a file
    G4int SetNtuple(G4TRNtupleDescription<TNTUPLE>* rntupleDescription);

    // Methods to bind ntuple (from base class)                
    using G4BaseRNtupleManager::SetNtupleIColumn; 
    using G4BaseRNtupleManager::SetNtupleFColumn; 
    using G4BaseRNtupleManager::SetNtupleDColumn; 
    using G4BaseRNtupleManager::SetNtupleSColumn; 

    // Methods to bind ntuple
    virtual G4bool SetNtupleIColumn(G4int ntupleId, 
                            const G4String& columnName, G4int& value) final;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, 
                            const G4String& columnName, G4float& value) final;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, 
                            const G4String& columnName, G4double& value) final;
    virtual G4bool SetNtupleSColumn(G4int ntupleId, 
                            const G4String& columnName, G4String& value) final;

    // Bind the ntuple columns of vector type
    virtual G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4int>& vector) override;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4float>& vector) override;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4double>& vector) override;

    using G4BaseRNtupleManager::GetNtupleRow;
    virtual G4bool GetNtupleRow(G4int ntupleId) final;
    
    // Access methods
    virtual G4int GetNofNtuples() const final;
  
    // Utility method
    G4TRNtupleDescription<TNTUPLE>*  GetNtupleDescriptionInFunction(G4int id, 
                                         G4String function,
                                         G4bool warn = true) const;

  private:
    // Fuctions which are specific to output type
    //
    virtual G4bool GetTNtupleRow(G4TRNtupleDescription<TNTUPLE>* rntupleDescription) = 0;

    // Common implementation
    //
    template <typename T> 
    G4bool SetNtupleTColumn(G4int ntupleId, const G4String& name, 
                            T& value);

    template <typename T> 
    G4bool SetNtupleTColumn(G4int ntupleId, const G4String& name, 
                            std::vector<T>& vector);

    // data members
    std::vector<G4TRNtupleDescription<TNTUPLE>*> fNtupleDescriptionVector;
};    

#include "G4TRNtupleManager.icc"

#endif

