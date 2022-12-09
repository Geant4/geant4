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

// The pure abstract base class for read Ntuple manager.
//
// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4VRNtupleManager_h
#define G4VRNtupleManager_h 1

#include "G4BaseAnalysisManager.hh"
#include "globals.hh"

#include <vector>

class G4VRNtupleManager : public G4BaseAnalysisManager
{
  // Disable using the object managers outside G4VAnalysisManager and
  // its messenger
  friend class G4VAnalysisReader;

  public:
    explicit G4VRNtupleManager(const G4AnalysisManagerState& state)
      : G4BaseAnalysisManager(state) {}
    G4VRNtupleManager() = delete;
    ~G4VRNtupleManager() override = default;

    // deleted copy constructor & assignment operator
    G4VRNtupleManager(const G4VRNtupleManager& rhs) = delete;
    G4VRNtupleManager& operator=(const G4VRNtupleManager& rhs) = delete;

  protected:
    // Methods to read ntuple from a file
    virtual G4int  ReadNtupleImpl(const G4String& ntupleName,
                            const G4String& fileName,
                            const G4String& dirName,
                            G4bool isUserFileName) = 0;

    // Methods for ntuple with id = FirstNtupleId
    virtual G4bool SetNtupleIColumn(const G4String& columnName,
                            G4int& value) = 0;
    virtual G4bool SetNtupleFColumn(const G4String& columnName,
                            G4float& value) = 0;
    virtual G4bool SetNtupleDColumn(const G4String& columnName,
                            G4double& value) = 0;
    virtual G4bool SetNtupleSColumn(const G4String& columnName,
                            G4String& value) = 0;
    // Bind the ntuple columns of vector type
    virtual G4bool SetNtupleIColumn(const G4String& columnName,
                            std::vector<G4int>& vector) = 0;
    virtual G4bool SetNtupleFColumn(const G4String& columnName,
                            std::vector<G4float>& vector) = 0;
    virtual G4bool SetNtupleDColumn(const G4String& columnName,
                            std::vector<G4double>& vector) = 0;
    virtual G4bool SetNtupleSColumn(const G4String& columnName,
                            std::vector<std::string>& vector) = 0;
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
    virtual G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName,
                            std::vector<G4int>& vector) = 0;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName,
                            std::vector<G4float>& vector) = 0;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName,
                            std::vector<G4double>& vector) = 0;
    virtual G4bool SetNtupleSColumn(G4int ntupleId, const G4String& columnName,
                            std::vector<std::string>& vector) = 0;
    virtual G4bool GetNtupleRow() = 0;
    virtual G4bool GetNtupleRow(G4int ntupleId) = 0;

    // Access methods
    virtual G4int GetNofNtuples() const = 0;
};

#endif

