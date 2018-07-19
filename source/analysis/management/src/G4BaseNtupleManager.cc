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

#include "G4BaseNtupleManager.hh"

//_____________________________________________________________________________
G4BaseNtupleManager::G4BaseNtupleManager(const G4AnalysisManagerState& state)
  : G4VNtupleManager(state),
    fFirstNtupleColumnId(0),
    fLockFirstNtupleColumnId(false)
{}

//_____________________________________________________________________________
G4BaseNtupleManager::~G4BaseNtupleManager()
{
}

//
// private methods
//

//_____________________________________________________________________________
G4int  G4BaseNtupleManager::GetCurrentNtupleId() const
{
  return  GetNofNtupleBookings() + fFirstId - 1;
}

//
// protected methods
//

//_____________________________________________________________________________
G4int G4BaseNtupleManager::CreateNtupleIColumn(const G4String& name, 
                                               std::vector<int>* vector)
{
  return CreateNtupleIColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4BaseNtupleManager::CreateNtupleFColumn(const G4String& name,
                                               std::vector<float>* vector)
{
  return CreateNtupleFColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4BaseNtupleManager::CreateNtupleDColumn(const G4String& name,
                                               std::vector<double>* vector)
{
  return CreateNtupleDColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4BaseNtupleManager::CreateNtupleSColumn(const G4String& name)
{
  return CreateNtupleSColumn(GetCurrentNtupleId(), name);
}

//_____________________________________________________________________________
void  G4BaseNtupleManager::FinishNtuple()
{
  return FinishNtuple(GetCurrentNtupleId());
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::FillNtupleIColumn(G4int id, G4int value)
{
  return FillNtupleIColumn(fFirstId, id, value);
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::FillNtupleFColumn(G4int id, G4float value)
{
  return FillNtupleFColumn(fFirstId, id, value);
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::FillNtupleDColumn(G4int id, G4double value)
{
  return FillNtupleDColumn(fFirstId, id, value);
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::FillNtupleSColumn(G4int id, const G4String& value)
{
  return FillNtupleSColumn(fFirstId, id, value);
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::AddNtupleRow()
{
  return AddNtupleRow(fFirstId);
}

//_____________________________________________________________________________
G4bool G4BaseNtupleManager::SetFirstNtupleColumnId(G4int firstId) 
{ 
  if ( fLockFirstNtupleColumnId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstNtupleColumnId as its value was already used.";
    G4Exception("G4BaseNtupleManager::SetFirstNtupleColumnId()",
                "Analysis_W013", JustWarning, description);
    return false;
  }              

  fFirstNtupleColumnId = firstId; 
  return true;
}

