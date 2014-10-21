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

// Manager class for Xml read ntuples.
// It implements functions specific to Xml read ntuples.
//
// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#ifndef G4XmlRNtupleManager_h
#define G4XmlRNtupleManager_h 1

#include "G4VRNtupleManager.hh"
#include "globals.hh"

#include "tools/raxml"

#include <vector>

struct G4XmlRNtupleDescription;

class G4XmlRNtupleManager : public G4VRNtupleManager
{
  friend class G4XmlAnalysisReader;

  protected:
    G4XmlRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4XmlRNtupleManager();

    // Methods to manipulate ntuples  
    G4bool IsEmpty() const;
    G4bool Reset();

    // Access methods
    tools::aida::ntuple* GetNtuple() const;
    tools::aida::ntuple* GetNtuple(G4int ntupleId) const;

    // Functions independent from the output type 
    //
    // Methods to read ntuple from a file
    G4int SetNtuple(G4XmlRNtupleDescription* rntupleDescription);
    // Methods for ntuple with id = FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(const G4String& columnName, 
                            G4int& value);
    virtual G4bool SetNtupleFColumn(const G4String& columnName, 
                            G4float& value);
    virtual G4bool SetNtupleDColumn(const G4String& columnName, 
                            G4double& value);
    virtual G4bool SetNtupleSColumn(const G4String& columnName, 
                            G4String& value);
    // Bind the ntuple columns of vector type
    virtual G4bool SetNtupleIColumn(const G4String& columnName, 
                            std::vector<G4int>& vector);
    virtual G4bool SetNtupleFColumn(const G4String& columnName, 
                            std::vector<G4float>& vector);
    virtual G4bool SetNtupleDColumn(const G4String& columnName, 
                            std::vector<G4double>& vector);
    // Methods for ntuple with id > FirstNtupleId                     
    virtual G4bool SetNtupleIColumn(G4int ntupleId, 
                            const G4String& columnName, G4int& value);
    virtual G4bool SetNtupleFColumn(G4int ntupleId, 
                            const G4String& columnName, G4float& value);
    virtual G4bool SetNtupleDColumn(G4int ntupleId, 
                            const G4String& columnName, G4double& value);
    virtual G4bool SetNtupleSColumn(G4int ntupleId, 
                            const G4String& columnName, G4String& value);
    // Bind the ntuple columns of vector type
    virtual G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4int>& vector);
    virtual G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4float>& vector);
    virtual G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4double>& vector);
    virtual G4bool GetNtupleRow();
    virtual G4bool GetNtupleRow(G4int ntupleId) ;
    
    // Access methods
    virtual G4int GetNofNtuples() const;
  
  private:
    // methods
    G4XmlRNtupleDescription*  GetNtupleInFunction(G4int id, 
                                         G4String function,
                                         G4bool warn = true) const;

    // data members
    std::vector<G4XmlRNtupleDescription*> fNtupleVector;
};    

// inline functions

inline G4int G4XmlRNtupleManager::GetNofNtuples() const
{ return fNtupleVector.size(); }

#endif

