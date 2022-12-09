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

// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#include "G4XmlRNtupleManager.hh"
#include "G4XmlRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;
using std::to_string;

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

//
// private methods
//

//_____________________________________________________________________________
G4int G4XmlRNtupleManager::ReadNtupleImpl(const G4String& ntupleName,
                                          const G4String& fileName,
                                          const G4String& /*dirName*/,
                                          G4bool isUserFileName)
{
  Message(kVL4, "read", "ntuple", ntupleName);

  // Ntuples are saved per object and per thread
  // but apply the ntuple name and the thread suffixes
  // only if fileName is not provided explicitly
  auto fullFileName = fileName;
  if ( ! isUserFileName ) {
    fullFileName = fFileManager->GetNtupleFileName(ntupleName);
  }

  auto handler = fFileManager->GetHandler<tools::aida::ntuple>(
                                 fullFileName, ntupleName, "ReadNtupleImpl");
  if (handler == nullptr) return kInvalidId;

  auto rntuple = static_cast<tools::aida::ntuple*>(handler->object());
  auto id = SetNtuple(new G4TRNtupleDescription<tools::aida::ntuple>(rntuple));

  Message(kVL2, "read", "ntuple", ntupleName, id > kInvalidId);

  return id;
}

//_____________________________________________________________________________
G4bool G4XmlRNtupleManager::SetNtupleIColumn(G4int ntupleId,
                                             const G4String& columnName,
                                             std::vector<G4int>& vector)
{
// Override base class default implementation

  Message(kVL4, "set", "ntuple I column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleIColumn");
  if (ntupleDescription == nullptr) return false;

  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fIVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

  Message(kVL2, "set", "ntuple I column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlRNtupleManager::SetNtupleFColumn(G4int ntupleId,
                                             const G4String& columnName,
                                             std::vector<G4float>& vector)
{
// Override base class default implementation

  Message(kVL4, "set", "ntuple F column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);


  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleFColumn");
  if (ntupleDescription == nullptr) return false;

  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fFVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

  Message(kVL4, "set", "ntuple F column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlRNtupleManager::SetNtupleDColumn(G4int ntupleId,
                                             const G4String& columnName,
                                             std::vector<G4double>& vector)
{
// Override base class default implementation

  Message(kVL4, "set", "ntuple D column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleDColumn");
  if (ntupleDescription == nullptr) return false;

  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fDVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

  Message(kVL2, "set", "ntuple D column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlRNtupleManager::SetNtupleSColumn(G4int ntupleId,
                                             const G4String& columnName,
                                             std::vector<std::string>& vector)
{
// Override base class default implementation

  Message(kVL4, "set", "ntuple S column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetNtupleSColumn");
  if (ntupleDescription == nullptr) return false;

  // not supported
  //tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  //ntupleBinding->add_column(columnName, vector);

  auto subNtuple = new tools::aida::ntuple(G4cout, columnName);
  ntupleDescription->fSVectorBindingMap[subNtuple] = &vector;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column_cid(columnName, *subNtuple);

  Message(kVL2, "set", "ntuple S column",
    " ntupleId " + to_string(ntupleId) + " " + columnName);

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
      Warn("Ntuple initialization failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  G4bool next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      Warn("Ntuple get_row() failed !!", fkClass, "GetTNtupleRow");
      return false;
    }

    // fill vector from sub ntuples
    for ( auto [key, value] : ntupleDescription->fIVectorBindingMap) {
      tools::aida::to_vector<int>(*key, *value);
    }
    for ( auto [key, value] : ntupleDescription->fFVectorBindingMap) {
      tools::aida::to_vector<float>(*key, *value);
    }
    for ( auto [key, value] : ntupleDescription->fDVectorBindingMap) {
      tools::aida::to_vector<double>(*key, *value);
    }
    for ( auto [key, value] : ntupleDescription->fSVectorBindingMap) {
      tools::aida::to_vector<std::string>(*key, *value);
    }
  }

  return next;
}
