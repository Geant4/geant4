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

#include "G4XmlRNtupleManager.hh"
#include "G4AnalysisManagerState.hh"

// 
// utility function (to be provided in tools)
//

namespace tools {
namespace aida {
template <class T>
bool to_vector(base_ntu& a_ntu,std::vector<T>& a_vec) {
  a_vec.clear();
  const std::vector<base_col*>& cols = a_ntu.cols();
  if(cols.empty()) return false;
  base_col* _base_col = cols.front();
  aida_col<T>* _col = safe_cast<base_col, aida_col<T> >(*_base_col);
  if(!_col) return false;
  a_ntu.start();
  uint64 _rows = a_ntu.rows();
  a_vec.resize(_rows);
  T v;
 {for(uint64 row=0;row<_rows;row++) {
    if(!a_ntu.next()) {a_vec.clear();return false;}
    if(!_col->get_entry(v)) {a_vec.clear();return false;}
    a_vec[row] = v;
  }}
  return true;
}
}}

//_____________________________________________________________________________
G4XmlRNtupleManager::G4XmlRNtupleManager(const G4AnalysisManagerState& state)
 : G4TRNtupleManager<tools::aida::ntuple>(state)
{}

//_____________________________________________________________________________
G4XmlRNtupleManager::~G4XmlRNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4XmlRNtupleManager::SetNtupleIColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4int>& vector)
{
// Override base class default implementation

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple I column", description);
  }  
#endif

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleIColumn");
  if ( ! ntupleDescription )  return false;   
  
  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fIVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

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
G4bool G4XmlRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4float>& vector)
{
// Override base class default implementation

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple F column of vector", description);
  }  
#endif

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleFColumn");
  if ( ! ntupleDescription )  return false;   
  
  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fFVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

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
G4bool G4XmlRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4double>& vector)
{
// Override base class default implementation

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple D column of vector", description);
  }  
#endif

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleDColumn");
  if ( ! ntupleDescription )  return false;   
  
  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fDVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

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
G4bool G4XmlRNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<tools::aida::ntuple>* ntupleDescription)
{
  auto ntuple = ntupleDescription->fNtuple;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {
    tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
    if ( ! ntuple->set_binding(std::cout, *ntupleBinding) ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Ntuple initialization failed !!"; 
      G4Exception("G4XmlRNtuple::GetNtupleRow()",
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
      G4Exception("G4XmlRNtuple::GetNtupleRow()",
                  "Analysis_WR021", JustWarning, description);
      return false;
    }

    // fill vector from sub ntuples

    {std::map<tools::aida::ntuple*, std::vector<int>* >::iterator it;
    for ( it = ntupleDescription->fIVectorBindingMap.begin(); 
          it != ntupleDescription->fIVectorBindingMap.end(); it++) {
      tools::aida::to_vector<int>(*(it->first), *(it->second));        
    }}
    {std::map<tools::aida::ntuple*, std::vector<float>* >::iterator it;
    for ( it = ntupleDescription->fFVectorBindingMap.begin(); 
          it != ntupleDescription->fFVectorBindingMap.end(); it++) {
      tools::aida::to_vector<float>(*(it->first), *(it->second));        
    }}
    {std::map<tools::aida::ntuple*, std::vector<double>* >::iterator it;
    for ( it = ntupleDescription->fDVectorBindingMap.begin(); 
          it != ntupleDescription->fDVectorBindingMap.end(); it++) {
      tools::aida::to_vector<double>(*(it->first), *(it->second));        
    }}
  }  

  return next;
}   
