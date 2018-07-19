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

// The base class for Ntuple managers. 
// It implements commen functions independent from the output type. 
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4BaseNtupleManager_h
#define G4BaseNtupleManager_h 1

#include "G4VNtupleManager.hh"
#include "globals.hh"

class G4BaseNtupleManager : public G4VNtupleManager
{
  public:
    explicit G4BaseNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4BaseNtupleManager();

    // deleted copy constructor & assignment operator
    G4BaseNtupleManager(const G4BaseNtupleManager& rhs) = delete;
    G4BaseNtupleManager& operator=(const G4BaseNtupleManager& rhs) =delete;

  protected:
    // Methods for handling ntuples
    virtual G4int CreateNtuple(const G4String& name, const G4String& title) = 0;

    // Create columns in the last created ntuple
    virtual G4int CreateNtupleIColumn(const G4String& name, 
                              std::vector<int>* vector) final;
    virtual G4int CreateNtupleFColumn(const G4String& name,
                              std::vector<float>* vector)  final;
    virtual G4int CreateNtupleDColumn(const G4String& name,
                              std::vector<double>* vector) final;
    virtual G4int CreateNtupleSColumn(const G4String& name);
    virtual void  FinishNtuple()  final;   

    // Create columns in the ntuple with given id
    virtual G4int CreateNtupleIColumn(G4int ntupleId, const G4String& name,
                                      std::vector<int>* vector) = 0;
    virtual G4int CreateNtupleFColumn(G4int ntupleId, const G4String& name,
                                      std::vector<float>* vector) = 0;
    virtual G4int CreateNtupleDColumn(G4int ntupleId, const G4String& name,
                                      std::vector<double>* vector) = 0;
    virtual G4int CreateNtupleSColumn(G4int ntupleId, const G4String& name) = 0;
    virtual void  FinishNtuple(G4int ntupleId) = 0; 
        
    // The ntuple column ids are generated automatically starting from 0; 
    // with the following function it is possible to change it 
    // to start from another value
    virtual G4bool SetFirstNtupleColumnId(G4int firstId) final; 
    G4int  GetFirstNtupleColumnId() const final;

    // Methods to fill ntuples
    // Methods for ntuple with id = FirstNtupleId                     
    virtual G4bool FillNtupleIColumn(G4int id, G4int value) final;
    virtual G4bool FillNtupleFColumn(G4int id, G4float value) final;
    virtual G4bool FillNtupleDColumn(G4int id, G4double value) final;
    virtual G4bool FillNtupleSColumn(G4int id, const G4String& value) final;
    virtual G4bool AddNtupleRow() final;

    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)                      
    virtual G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value) = 0;
    virtual G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value) = 0;
    virtual G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value) = 0;
    virtual G4bool FillNtupleSColumn(G4int ntupleId, G4int columnId, 
                                     const G4String& value) = 0;
    virtual G4bool AddNtupleRow(G4int ntupleId) = 0;
   
  protected:
    G4int   fFirstNtupleColumnId;
    G4bool  fLockFirstNtupleColumnId;

  private:
    // methods
    G4int GetCurrentNtupleId() const;
};
// inline functions

inline G4int G4BaseNtupleManager::GetFirstNtupleColumnId() const {
  return fFirstNtupleColumnId;
}  
    
#endif

