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
// $Id: G4RootRNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#include "G4RootRNtupleManager.hh"
#include "G4RootRNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"

//_____________________________________________________________________________
G4RootRNtupleManager::G4RootRNtupleManager(const G4AnalysisManagerState& state)
 : G4VRNtupleManager(state),
   fNtupleVector()
{
}

//_____________________________________________________________________________
G4RootRNtupleManager::~G4RootRNtupleManager()
{  
  std::vector<G4RootRNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    delete (*it);
  }   
}

// 
// private methods
//

//_____________________________________________________________________________
G4RootRNtupleDescription* G4RootRNtupleManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootRNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_WR011", JustWarning, description);
    }
    return nullptr;         
  }
  
  return fNtupleVector[index];
}
  
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::IsEmpty() const
{
  return ! fNtupleVector.size();
}  
 
//_____________________________________________________________________________
G4bool G4RootRNtupleManager::Reset()
{
// Reset ntuples

  std::vector<G4RootRNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    // ntuple is deleted automatically when file is closed
    // delete (*it)->fNtuple;
    (*it)->fNtuple=0; 
  }  
  
  return true;
}  

//_____________________________________________________________________________
tools::rroot::ntuple* G4RootRNtupleManager::GetNtuple() const
{
  return GetNtuple(fFirstId);
}  

//_____________________________________________________________________________
tools::rroot::ntuple* G4RootRNtupleManager::GetNtuple(G4int ntupleId) const
{
  G4RootRNtupleDescription* rntupleDescription
    = GetNtupleInFunction(ntupleId, "GetRNtuple");

  if ( ! rntupleDescription ) return nullptr; 
    
  return rntupleDescription->fNtuple;  
}  

//_____________________________________________________________________________
G4int G4RootRNtupleManager::SetNtuple(G4RootRNtupleDescription* rntupleDescription)
{
  G4int id = fNtupleVector.size() + fFirstId;

  fNtupleVector.push_back(rntupleDescription);
  
  return id;
}  

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                             G4int& value)
{
  return SetNtupleIColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                             G4float& value)
{                                             
  return SetNtupleFColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                             G4double& value)
{                                             
  return SetNtupleDColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleSColumn(const G4String& columnName, 
                                             G4String& value)
{                                             
  return SetNtupleSColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                             std::vector<G4int>& vector)
{
  return SetNtupleIColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                             std::vector<G4float>& vector)
{                                             
  return SetNtupleFColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                             std::vector<G4double>& vector)
{                                             
  return SetNtupleDColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleIColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             G4int& value)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple I column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleIColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, value);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) { 
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple I colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             G4float& value)
{                                             
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple F column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleFColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, value);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple F colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             G4double& value)
{                                             
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple D column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleDColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, value);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple D colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleSColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             G4String& value)
{                                             
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple S column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleSColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, value);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple S colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleIColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4int>& vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple I column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleIColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple I colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4float>& vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple F column of vector", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleFColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple F colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4double>& vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple D column of vector", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleDColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple D colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetNtupleRow()
{ 
  return GetNtupleRow(fFirstId);                                            
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetNtupleRow(G4int ntupleId)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("get", "ntuple row", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtupleRow");
  if ( ! ntupleDescription )  return false;   
  
  tools::rroot::ntuple* ntuple = ntupleDescription->fNtuple;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      G4ExceptionDescription description;
      description << " ntupleId " << ntupleId;  
      fState.GetVerboseL4()->Message("initialize", "ntuple", description);
    }  
#endif
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple initialization failed !!"; 
      G4Exception("G4RootRNtuple::GetNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  G4bool next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple get_row() failed !!"; 
      G4Exception("G4RootRNtuple::GetNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }
  }  

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL2()->Message("get", "ntuple row", description, true);
  }  
#endif

  return next;
}   
