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

#include "G4TRNtupleManager.hh"
#include "globals.hh"

#include "tools/raxml"

#include <vector>

struct G4XmlRNtupleDescription;

class G4XmlRNtupleManager : public G4TRNtupleManager<tools::aida::ntuple>
{
  friend class G4XmlAnalysisReader;

  protected:
    G4XmlRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4XmlRNtupleManager();

    // Methods from base classes            
    using G4BaseRNtupleManager::SetNtupleIColumn; 
    using G4BaseRNtupleManager::SetNtupleFColumn; 
    using G4BaseRNtupleManager::SetNtupleDColumn; 
    using G4TRNtupleManager<tools::aida::ntuple>::SetNtupleIColumn; 
    using G4TRNtupleManager<tools::aida::ntuple>::SetNtupleFColumn; 
    using G4TRNtupleManager<tools::aida::ntuple>::SetNtupleDColumn; 

    // Override base class functions for vector columns
    // in order to fill the maps
  
    virtual G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4int>& vector) final;
    virtual G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4float>& vector) final;
    virtual G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<G4double>& vector) final;

  private:
    // Methods from the templated base class
    //
    virtual G4bool GetTNtupleRow(G4TRNtupleDescription<tools::aida::ntuple>* ntupleDescription) final;
};    

#endif

