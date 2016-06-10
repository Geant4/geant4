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

// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#include "G4CsvRNtupleManager.hh"
#include "G4CsvRNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"

//_____________________________________________________________________________
G4CsvRNtupleManager::G4CsvRNtupleManager(const G4AnalysisManagerState& state)
 : G4VRNtupleManager(state),
   fNtupleVector()
{
}

//_____________________________________________________________________________
G4CsvRNtupleManager::~G4CsvRNtupleManager()
{  
  std::vector<G4CsvRNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    delete (*it);
  }   
}

// 
// private methods
//

//_____________________________________________________________________________
G4CsvRNtupleDescription* G4CsvRNtupleManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4CsvRNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_WR011", JustWarning, description);
    }
    return 0;         
  }
  
  return fNtupleVector[index];
}
  
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::IsEmpty() const
{
  return ! fNtupleVector.size();
}  
 
//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::Reset()
{
// Reset ntuples

  std::vector<G4CsvRNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    // ntuple is deleted automatically when file is closed
    // delete (*it)->fNtuple;
    (*it)->fNtuple=0; 
  }  
  
  return true;
}  

//_____________________________________________________________________________
tools::rcsv::ntuple* G4CsvRNtupleManager::GetNtuple() const
{
  return GetNtuple(fFirstId);
}  

//_____________________________________________________________________________
tools::rcsv::ntuple* G4CsvRNtupleManager::GetNtuple(G4int ntupleId) const
{
  G4CsvRNtupleDescription* rntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");

  if ( ! rntupleDescription ) return 0; 

  return rntupleDescription->fNtuple;  
}  

//_____________________________________________________________________________
G4int G4CsvRNtupleManager::SetNtuple(G4CsvRNtupleDescription* rntupleDescription)
{
  G4int id = fNtupleVector.size() + fFirstId;

  fNtupleVector.push_back(rntupleDescription);
  
  return id;
}  

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                             G4int& value)
{
  return SetNtupleIColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                             G4float& value)
{                                             
  return SetNtupleFColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                             G4double& value)
{                                             
  return SetNtupleDColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleIColumn(const G4String& columnName, 
                                             std::vector<G4int>& vector)
{
  return SetNtupleIColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleFColumn(const G4String& columnName, 
                                             std::vector<G4float>& vector)
{                                             
  return SetNtupleFColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleDColumn(const G4String& columnName, 
                                             std::vector<G4double>& vector)
{                                             
  return SetNtupleDColumn(fFirstId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleSColumn(const G4String& columnName, 
                                             G4String& value)
{                                             
  return SetNtupleSColumn(fFirstId, columnName, value);
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::SetNtupleIColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleIColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::SetNtupleSColumn(G4int ntupleId, 
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

  G4CsvRNtupleDescription* ntupleDescription
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
G4bool G4CsvRNtupleManager::GetNtupleRow()
{ 
  return GetNtupleRow(fFirstId);                                            
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::GetNtupleRow(G4int ntupleId)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("get", "ntuple row", description);
  }  
#endif

  G4CsvRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtupleRow");
  if ( ! ntupleDescription )  return false;   
  
  tools::rcsv::ntuple* ntuple = ntupleDescription->fNtuple;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {
    tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple initialization failed !!"; 
      G4Exception("G4CsvRNtuple::GetNtupleRow()",
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
      G4Exception("G4CsvRNtuple::GetNtupleRow()",
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
