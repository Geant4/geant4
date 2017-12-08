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

// The base class for read Ntuple manager. 
// It defines functions independent from the output type. 
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4BaseRNtupleManager_h
#define G4BaseRNtupleManager_h 1

#include "G4VRNtupleManager.hh"
#include "globals.hh"

#include <vector>

class G4BaseRNtupleManager : public G4VRNtupleManager
{
  // Disable using the object managers outside G4VAnalysisManager and
  // its messenger
  friend class G4VAnalysisReader;

  public:
    explicit G4BaseRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4BaseRNtupleManager();

    // deleted copy constructor & assignment operator
    G4BaseRNtupleManager(const G4BaseRNtupleManager& rhs) = delete;
    G4BaseRNtupleManager& operator=(const G4BaseRNtupleManager& rhs) = delete;

  protected:   
    // Methods to read ntuple from a file
    // Methods for ntuple with id = FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(const G4String& columnName, 
                            G4int& value) final;
    virtual G4bool SetNtupleFColumn(const G4String& columnName, 
                            G4float& value) final;
    virtual G4bool SetNtupleDColumn(const G4String& columnName, 
                            G4double& value) final;
    virtual G4bool SetNtupleSColumn(const G4String& columnName, 
                            G4String& value) final;
    // Methods for ntuple with id > FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(G4int ntupleId, 
                            const G4String& columnName, G4int& value)= 0;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, 
                            const G4String& columnName, G4float& value)= 0;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, 
                            const G4String& columnName, G4double& value)= 0;
    virtual G4bool SetNtupleSColumn(G4int ntupleId, 
                            const G4String& columnName, G4String& value)= 0;
    // Bind the ntuple columns of vector type
    // Methods for ntuple with id = FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(const G4String& columnName, 
                            std::vector<G4int>& vector) final;
    virtual G4bool SetNtupleFColumn(const G4String& columnName, 
                            std::vector<G4float>& vector) final;
    virtual G4bool SetNtupleDColumn(const G4String& columnName, 
                            std::vector<G4double>& vector) final;
    // Methods for ntuple with id > FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4int>& vector) = 0;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4float>& vector) = 0;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4double>& vector) = 0;
    virtual G4bool GetNtupleRow() final;
    virtual G4bool GetNtupleRow(G4int ntupleId) = 0;
    
    // Access methods
    virtual G4int GetNofNtuples() const = 0;

  private:
    // methods
    G4int GetCurrentNtupleId() const;
};
    
#endif

