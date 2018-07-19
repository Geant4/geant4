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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4BaseRNtupleManager.hh"

//_____________________________________________________________________________
G4BaseRNtupleManager::G4BaseRNtupleManager(const G4AnalysisManagerState& state)
  : G4VRNtupleManager(state)
{
}

//_____________________________________________________________________________
G4BaseRNtupleManager::~G4BaseRNtupleManager()
{
}

//
// private methods
//
//_____________________________________________________________________________
G4int  G4BaseRNtupleManager::GetCurrentNtupleId() const
{
  return  GetNofNtuples() + fFirstId - 1;
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                           G4int& value)
{
  return SetNtupleIColumn(GetCurrentNtupleId(), columnName, value);
}
//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                           G4float& value)
{
  return SetNtupleFColumn(GetCurrentNtupleId(), columnName, value);  
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                           G4double& value)
{
  return SetNtupleDColumn(GetCurrentNtupleId(), columnName, value);    
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleSColumn(const G4String& columnName, 
                                           G4String& value)
{
  return SetNtupleSColumn(GetCurrentNtupleId(), columnName, value);      
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                           std::vector<G4int>& vector)
{
  return SetNtupleIColumn(GetCurrentNtupleId(), columnName, vector);
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                           std::vector<G4float>& vector)
{
  return SetNtupleFColumn(GetCurrentNtupleId(), columnName, vector);
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                           std::vector<G4double>& vector)
{
  return SetNtupleDColumn(GetCurrentNtupleId(), columnName, vector);
}

//_____________________________________________________________________________
G4bool G4BaseRNtupleManager::GetNtupleRow()
{
  return GetNtupleRow(fFirstId);  
}
