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
//
/// \file hbook/src/ExG4HbookNtupleManager.cc
/// \brief Implementation of the ExG4HbookNtupleManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookNtupleManager.hh"
#include "ExG4HbookFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4UnitsTable.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
ExG4HbookNtupleManager::ExG4HbookNtupleManager(const G4AnalysisManagerState& state)
 : G4VNtupleManager(state),
   fNtupleHbookIdOffset(-1),
   fNtupleDescriptionVector(),
   fNtupleVector()
{
}

//_____________________________________________________________________________
ExG4HbookNtupleManager::~ExG4HbookNtupleManager()
{  
  // Reset();

  std::vector<ExG4HbookNtupleDescription*>::iterator it;  
  for (it = fNtupleDescriptionVector.begin(); it != fNtupleDescriptionVector.end(); it++ ) {
    delete *it;
  }   
}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookNtupleManager::SetNtupleHbookIdOffset()
{
// Set  fH1HbookIdOffset if needed

  if ( fNtupleHbookIdOffset == -1 ) {
    if ( fFirstId > 0 ) 
      fNtupleHbookIdOffset = 0;
    else
      fNtupleHbookIdOffset = 1;
        
    if ( fNtupleHbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "Ntuple will be defined in HBOOK with ID = G4_firstNtupleId + 1";
      G4Exception("ExG4HbookNtupleManager::SetNtupleHbookIdOffset()",
                  "Analysis_W013", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( ! fNtupleDescriptionVector.size() ) return;     
  
  // Set fNtupleHbookIdOffset if needed
  SetNtupleHbookIdOffset();
  
  G4int index = 0;
  std::vector<ExG4HbookNtupleDescription*>::iterator itn;  
  for (itn = fNtupleDescriptionVector.begin(); itn != fNtupleDescriptionVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;
    if ( ! ntupleBooking ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "ntuple", ntupleBooking->name());
#endif

    // Create an "ntuple" directory both in memory and in the file
    fFileManager->CreateNtupleDirectory();
    G4int hbookIndex = fNtupleHbookIdOffset + index + fFirstId;
    ++index;

    // We should be under //PAWC/LUN1/ntuple
    (*itn)->fNtuple
      = new tools::hbook::wntuple(hbookIndex, G4cout, *ntupleBooking);
    fNtupleVector.push_back((*itn)->fNtuple);  

    if ( ntupleBooking->columns().size() ) {
      // store ntuple columns in local maps
      const std::vector<tools::column_booking>& columns 
        = ntupleBooking->columns();
      std::vector<tools::column_booking>::const_iterator it;
      G4int counter = 0;
      for ( it = columns.begin(); it!=columns.end(); ++it) {
        if ( it->cls_id() == tools::_cid(int(0) ) ) {
          (*itn)->fNtupleIColumnMap[counter++] 
            = (*itn)->fNtuple->find_column<int>(it->name());
        }
        else if( it->cls_id() == tools::_cid(float(0) ) ) {
          (*itn)->fNtupleFColumnMap[counter++] 
            = (*itn)->fNtuple->find_column<float>(it->name());
        } 
        else if(it->cls_id()== tools::_cid(double(0))) {
          (*itn)->fNtupleDColumnMap[counter++] 
            = (*itn)->fNtuple->find_column<double>(it->name());
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << it->name();
          G4Exception("G4HbookAnalysisManager::CreateNtupleFromBooking()",
                      "Analysis_W002", JustWarning, description);
        }
      }
    }
    FinishNtuple();
#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()
      ->Message("create from booking", "ntuple", ntupleBooking->name());
#endif
  }  
}   

//_____________________________________________________________________________
tools::hbook::wntuple::column<int>*    
ExG4HbookNtupleManager::GetNtupleIColumn(G4int ntupleId, G4int columnId) const
{
  ExG4HbookNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::hbook::wntuple::column<int>* >& ntupleIColumnMap
    = ntupleDecription->fNtupleIColumnMap;
  std::map<G4int, tools::hbook::wntuple::column<int>* >::const_iterator it
    = ntupleIColumnMap.find(columnId);
  if ( it == ntupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleIColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::hbook::wntuple::column<float>*  
ExG4HbookNtupleManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  ExG4HbookNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleFColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::hbook::wntuple::column<float>* >& ntupleFColumnMap
    = ntupleDecription->fNtupleFColumnMap;
  std::map<G4int, tools::hbook::wntuple::column<float>* >::const_iterator it
    = ntupleFColumnMap.find(columnId);
  if ( it == ntupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleFColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  

//_____________________________________________________________________________
tools::hbook::wntuple::column<double>* 
ExG4HbookNtupleManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  ExG4HbookNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleDColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::hbook::wntuple::column<double>* >& ntupleDColumnMap
    = ntupleDecription->fNtupleDColumnMap;
  std::map<G4int, tools::hbook::wntuple::column<double>* >::const_iterator it
    = ntupleDColumnMap.find(columnId);
  if ( it == ntupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleDColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
void ExG4HbookNtupleManager::Reset()
{
// Reset ntuple  

  std::vector<ExG4HbookNtupleDescription*>::iterator it3;  
  for (it3 = fNtupleDescriptionVector.begin(); it3 != fNtupleDescriptionVector.end(); it3++ ) {
    delete (*it3)->fNtuple;
    (*it3)->fNtuple = 0;
  }  
  fNtupleVector.clear(); 
}  
 
//
// protected methods
//

//_____________________________________________________________________________
ExG4HbookNtupleDescription* ExG4HbookNtupleManager::GetNtupleInFunction(
                                      G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool /*onlyIfActive*/) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4HbookAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }
  
  return fNtupleDescriptionVector[index];
}
  
// 
// public methods
//

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtuple(const G4String& name, 
                                           const G4String& title)
{
  // Create an "ntuple" directory both in memory and in the file
  if ( fFileManager->IsFile() ) 
    fFileManager->CreateNtupleDirectory();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple", name);
#endif

  // Create ntuple description
  G4int index = fNtupleDescriptionVector.size();
  ExG4HbookNtupleDescription* ntupleDescription
    = new ExG4HbookNtupleDescription();
  fNtupleDescriptionVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking 
    = new tools::ntuple_booking(name, title);
           // ntuple booking object is deleted in destructor

  // Set fNtupleHbookIdOffset if needed
  SetNtupleHbookIdOffset();
  
  // Create ntuple if the file is open
  // We should be under //PAWC/LUN1/ntuple
  if ( fFileManager->IsFile() ) {
    G4int hbookIndex = fNtupleHbookIdOffset + index + fFirstId;
    ntupleDescription->fNtuple 
      = new tools::hbook::wntuple(hbookIndex, name);
           // ntuple object is deleted when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);       
  }  

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << index + fFirstId; 
    fState.GetVerboseL2()->Message("create", "ntuple", description);
  }  
#endif

  return index + fFirstId;
}                                         

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleIColumn(const G4String& name,
                                                  std::vector<int>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleIColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleFColumn(const G4String& name,
                                                  std::vector<float>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleFColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleDColumn(const G4String& name,
                                                  std::vector<double>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleDColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleSColumn(const G4String& name)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleSColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void ExG4HbookNtupleManager::FinishNtuple()
{ 
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  FinishNtuple(ntupleId);
}
  
//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleIColumn(G4int ntupleId, 
                                                    const G4String& name,
                                                    std::vector<int>* vector)
{
  // Vector columns are not supported in HBOOK
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
      << "Vector columns are not supported in HBOOK."; 
    G4Exception("(ExG4HbookNtupleManager::CreateNtupleIColumn)",
                "Analysis_W002", FatalException, description);
    return kInvalidId;            
  }
                
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple I column", description);
  }  
#endif

  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleIColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleIColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::hbook::wntuple::column<int>* column 
      = ntupleDescription->fNtuple->create_column<int>(name);  
    ntupleDescription->fNtupleIColumnMap[index] = column;
  }  
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL2()->Message("create", "ntuple I column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleFColumn(G4int ntupleId, 
                                                  const G4String& name,
                                                  std::vector<float>* vector)
{ 
  // Vector columns are not supported in HBOOK
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
      << "Vector columns are not supported in HBOOK."; 
    G4Exception("(ExG4HbookNtupleManager::CreateNtupleFColumn)",
                "Analysis_W002", FatalException, description);
    return kInvalidId;            
  }
                
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple F column", description);
  } 
#endif

  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleFColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleFColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::hbook::wntuple::column<float>* column 
      = ntupleDescription->fNtuple->create_column<float>(name);  
    ntupleDescription->fNtupleFColumnMap[index] = column;
  }
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL2()->Message("create", "ntuple F column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         


//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleDColumn(G4int ntupleId, 
                                                  const G4String& name,
                                                  std::vector<double>* vector)
{
  // Vector columns are not supported in HBOOK
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
      << "Vector columns are not supported in HBOOK."; 
    G4Exception("(ExG4HbookNtupleManager::CreateNtupleDColumn)",
                "Analysis_W002", FatalException, description);
    return kInvalidId;            
  }
                
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple D column", description);
  }  
#endif

  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::hbook::wntuple::column<double>* column 
      = ntupleDescription->fNtuple->create_column<double>(name);  
    ntupleDescription->fNtupleDColumnMap[index] = column;
  }
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL2()->Message("create", "ntuple D column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int ExG4HbookNtupleManager::CreateNtupleSColumn(G4int /*ntupleId*/, 
                                                  const G4String& /*name*/)
{
  G4ExceptionDescription description;
  description << "      " 
    << "Columns of string type are not supported in HBOOK."; 
  G4Exception("(ExG4HbookNtupleManager::CreateNtupleDColumn)",
              "Analysis_W002", FatalException, description);
  return kInvalidId;
}

//_____________________________________________________________________________
void ExG4HbookNtupleManager::FinishNtuple(G4int ntupleId)
{ 
  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  tools::hbook::wntuple* ntuple = ntupleDescription->fNtuple;  

  if ( ! ntuple ) return;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("finish", "ntuple", ntupleDescription->fNtupleBooking->name());
#endif

  // Return to //PAWC/LUN1 :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");

  //fNtuple->add_row_beg();
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()
      ->Message("finish", "ntuple", ntupleDescription->fNtupleBooking->name());
#endif
}
  
//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleIColumn(G4int columnId, G4int value)
{
  return FillNtupleIColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleFColumn(G4int columnId, G4float value)
{
  return FillNtupleFColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleDColumn(G4int columnId, G4double value)
{
  return FillNtupleDColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleSColumn(G4int columnId, 
                                                 const G4String& value)
{
  return FillNtupleSColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::AddNtupleRow()
{ 
  return AddNtupleRow(fFirstId);
}  

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleIColumn(
                                            G4int ntupleId, G4int columnId, 
                                            G4int value)
{
  tools::hbook::wntuple::column<int>* column 
    = GetNtupleIColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleIColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
 #ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "ntuple I column", description);
  }  
#endif
 return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleFColumn(
                                            G4int ntupleId, G4int columnId, 
                                            G4float value)
{
  tools::hbook::wntuple::column<float>* column 
    = GetNtupleFColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleFColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleDColumn(
                                            G4int ntupleId, G4int columnId, 
                                            G4double value)
{
  tools::hbook::wntuple::column<double>* column 
    = GetNtupleDColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleDColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::FillNtupleSColumn(
                                   G4int /*ntupleId*/, G4int /*columnId*/, 
                                   const G4String& /*value*/)
{
  G4ExceptionDescription description;
  description << "      " 
    << "Columns of string type are not supported in HBOOK."; 
  G4Exception("(ExG4HbookNtupleManager::FillNtupleSColumn)",
              "Analysis_W011", JustWarning, description);
  return kInvalidId;
}  

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::AddNtupleRow(G4int ntupleId)
{ 
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description);
  }  
#endif

  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "AddNtupleRow");

  if ( ! ntupleDescription || ! ntupleDescription->fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << " does not exist. ";
    G4Exception("G4HbookAnalysisManager::AddNtupleRow()",
                "Analysis_W011", JustWarning, description);
    return false;
  }  
  
  ntupleDescription->fNtuple->add_row();
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description, true);
  }  
#endif
  return true;
}
 
//_____________________________________________________________________________
tools::hbook::wntuple* ExG4HbookNtupleManager::GetNtuple() const
{
  return GetNtuple(fFirstId);
}

//_____________________________________________________________________________
tools::hbook::wntuple* ExG4HbookNtupleManager::GetNtuple(G4int ntupleId) const
{
  ExG4HbookNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");
    
  return ntupleDescription->fNtuple;  
}

//_____________________________________________________________________________
G4bool ExG4HbookNtupleManager::SetNtupleHbookIdOffset(G4int offset) 
{
  if ( fNtupleDescriptionVector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set NtupleHbookIdOffset as some ntuples already exist.";
    G4Exception("G4HbookAnalysisManager::SetNtupleHbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }
  
  if ( fFirstId + offset < 1 ) {
    G4ExceptionDescription description;
    description << "The first ntuple HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetNtupleHbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }
  
  fNtupleHbookIdOffset = offset;
  return true;
}  

#endif
