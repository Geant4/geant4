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
// $Id: G4CsvNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvNtupleManager.hh"
#include "G4CsvNtupleDescription.hh"
#include "G4CsvFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4CsvNtupleManager::G4CsvNtupleManager(const G4AnalysisManagerState& state)
 : G4VNtupleManager(state),
   fFileManager(0),
   fNtupleDescriptionVector(),
   fNtupleVector(),
   fIsCommentedHeader(true),
   fIsHippoHeader(false)
{
}

//_____________________________________________________________________________
G4CsvNtupleManager::~G4CsvNtupleManager()
{  
  std::vector<G4CsvNtupleDescription*>::iterator it;  
  for (it = fNtupleDescriptionVector.begin(); it != fNtupleDescriptionVector.end(); it++ ) {
    delete (*it);
  }   
}

// 
// private methods
//

//_____________________________________________________________________________
tools::wcsv::ntuple::column<int>*    
G4CsvNtupleManager::GetNtupleIColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wcsv::ntuple::column<int>* >& ntupleIColumnMap
    = ntupleDecription->fNtupleIColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<int>* >::const_iterator it
    = ntupleIColumnMap.find(columnId);
  if ( it == ntupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::GetNtupleIColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wcsv::ntuple::column<float>*  
G4CsvNtupleManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleFColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wcsv::ntuple::column<float>* >& ntupleFColumnMap
    = ntupleDecription->fNtupleFColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<float>* >::const_iterator it
    = ntupleFColumnMap.find(columnId);
  if ( it == ntupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::GetNtupleFColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wcsv::ntuple::column<double>* 
G4CsvNtupleManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleDColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wcsv::ntuple::column<double>* >& ntupleDColumnMap
    = ntupleDecription->fNtupleDColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<double>* >::const_iterator it
    = ntupleDColumnMap.find(columnId);
  if ( it == ntupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::GetNtupleDColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
tools::wcsv::ntuple::column<std::string>* 
G4CsvNtupleManager::GetNtupleSColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleSColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wcsv::ntuple::column<std::string>* >& ntupleSColumnMap
    = ntupleDecription->fNtupleSColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<std::string>* >::const_iterator it
    = ntupleSColumnMap.find(columnId);
  if ( it == ntupleSColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::GetNtupleSColumn()",
                "Analysis_W011", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
G4CsvNtupleDescription* G4CsvNtupleManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool /*onlyIfActive*/) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4CsvNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }
  
  return fNtupleDescriptionVector[index];
}
  
//_____________________________________________________________________________
G4bool G4CsvNtupleManager::WriteHeader(tools::wcsv::ntuple* ntuple) const
{
// Write header if ntuple already exists and if this option is activated.
// When both Hippo and Commented headers are selected, only Commented
// header, which reading is supported.
// Return false only if an error occurred. 

  if ( fIsCommentedHeader ) {
    return ntuple->write_commented_header(G4cout);
  }
  
  // write hippo header (if activated and if not commented header)
  if ( fIsHippoHeader ) {
    ntuple->write_hippo_header();
    return true;
  }
  
  return true;
}

// 
// protected methods
//

//_____________________________________________________________________________
void G4CsvNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  std::vector<G4CsvNtupleDescription*>::iterator itn;  
  for (itn = fNtupleDescriptionVector.begin(); itn != fNtupleDescriptionVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;  
    if ( ! ntupleBooking ) continue;
    
    // Do not create ntuple if it already exists
    if ( (*itn)->fNtuple ) continue;
    
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "ntuple", ntupleBooking->name());
#endif

    // create a file for this ntuple
    if ( ! fFileManager->CreateNtupleFile((*itn)) ) continue;

    // create ntuple
    (*itn)->fNtuple
      = new tools::wcsv::ntuple(*((*itn)->fFile), G4cerr, *ntupleBooking);
    fNtupleVector.push_back((*itn)->fNtuple);  
  
    // write header (if activated)
    if ( ! WriteHeader((*itn)->fNtuple) ) {
       G4ExceptionDescription description;
       description << "      " 
                   << "Writing ntuple header has failed. ";
       G4Exception("G4CsvNtupleManager::CreateNtupleFromBooking()",
                   "Analysis_W021", JustWarning, description);
    }

    if ( ntupleBooking->columns().size() ) {
      // store ntuple columns in local maps
      const std::vector<tools::column_booking>& columns 
        = ntupleBooking->columns();
      std::vector<tools::column_booking>::const_iterator it;
      G4int index = 0;
      for ( it = columns.begin(); it!=columns.end(); ++it) {
        if ( it->cls_id() == tools::_cid(int(0) ) ) {
          (*itn)->fNtupleIColumnMap[index++] 
            = (*itn)->fNtuple->find_column<int>(it->name());
        }
        else if ( it->cls_id() == tools::_cid(float(0) ) ) {
          (*itn)->fNtupleFColumnMap[index++] 
            = (*itn)->fNtuple->find_column<float>(it->name());
        } 
        else if ( it->cls_id() == tools::_cid(double(0))) {
          (*itn)->fNtupleDColumnMap[index++] 
            = (*itn)->fNtuple->find_column<double>(it->name());
        }
        else if ( it->cls_id() == tools::_cid(std::string(""))) {
          (*itn)->fNtupleSColumnMap[index++] 
            = (*itn)->fNtuple->find_column<std::string>(it->name());
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << it->cls_id();
          G4Exception("G4CsvNtupleManager::CreateNtupleFromBooking()",
                      "Analysis_W002", JustWarning, description);
        }
      }
    }
#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()
        ->Message("create from booking", "ntuple", ntupleBooking->name());
#endif
  }
}   

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::IsEmpty() const
{
  return ! fNtupleDescriptionVector.size();
}  
 
//_____________________________________________________________________________
G4bool G4CsvNtupleManager::Reset()
{
  std::vector<G4CsvNtupleDescription*>::iterator it;  
  for (it = fNtupleDescriptionVector.begin(); it != fNtupleDescriptionVector.end(); it++ ) {
    delete (*it)->fNtuple; 
    (*it)->fNtuple = 0;
  }  
  fNtupleVector.clear(); 
  
  return true;
}  
 
//_____________________________________________________________________________
tools::wcsv::ntuple* G4CsvNtupleManager::GetNtuple() const
{
  return GetNtuple(fFirstId);
}  

//_____________________________________________________________________________
tools::wcsv::ntuple* G4CsvNtupleManager::GetNtuple(G4int ntupleId) const
{
  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");

  if ( ! ntupleDescription ) return 0; 
    
  return ntupleDescription->fNtuple;  
}  

//_____________________________________________________________________________
G4int G4CsvNtupleManager::CreateNtuple(const G4String& name, 
                                       const G4String& title)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple", name);
#endif

  // Create ntuple description
  G4int index = fNtupleDescriptionVector.size();
  G4CsvNtupleDescription* ntupleDescription
    = new G4CsvNtupleDescription();
  fNtupleDescriptionVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking 
    = new tools::ntuple_booking(name, title);
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open (what means here that
  // a filename was already set)
  if ( fFileManager->GetFileName().size() ) {
    if ( fFileManager->CreateNtupleFile(ntupleDescription) ) {
      ntupleDescription->fNtuple 
        = new tools::wcsv::ntuple(*(ntupleDescription->fFile));
           // ntuple object is deleted when closing a file
      (ntupleDescription->fNtuple)->set_title(title); 
      fNtupleVector.push_back(ntupleDescription->fNtuple);       
    }       
  }  

  fLockFirstId = true;

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
G4int G4CsvNtupleManager::CreateNtupleIColumn(const G4String& name,
                                              std::vector<int>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleIColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4CsvNtupleManager::CreateNtupleFColumn(const G4String& name,
                                              std::vector<float>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleFColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4CsvNtupleManager::CreateNtupleDColumn(const G4String& name,
                                              std::vector<double>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleDColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4CsvNtupleManager::CreateNtupleSColumn(const G4String& name)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleSColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void G4CsvNtupleManager::FinishNtuple()
{ 
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  FinishNtuple(ntupleId);
}
   
//_____________________________________________________________________________
G4int G4CsvNtupleManager::CreateNtupleIColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<int>* vector)
{
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple columns of vector are not supported." ;
    G4Exception("G4CsvAnalysisManager::CreateNtupleIColumn", 
                "Analysis_W002", JustWarning, description);
    return 0;
  }                

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple I column", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleIColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvNtupleManager::CreateNtupleIColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<int>* column 
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
G4int G4CsvNtupleManager::CreateNtupleFColumn(G4int ntupleId, 
                                              const G4String& name,
                                              std::vector<float>* vector)
{
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple columns of vector are not supported." ;
    G4Exception("G4CsvAnalysisManager::CreateNtupleFColumn", 
                "Analysis_W002", JustWarning, description);
    return 0;
  }                

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple F column", description);
  } 
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleFColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvNtupleManager::CreateNtupleFColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<float>* column 
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
G4int G4CsvNtupleManager::CreateNtupleDColumn(G4int ntupleId,
                                              const G4String& name,
                                              std::vector<double>* vector)
{
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple columns of vector are not supported." ;
    G4Exception("G4CsvAnalysisManager::CreateNtupleDColumn", 
                "Analysis_W002", JustWarning, description);
    return 0;
  }                

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple D column", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvNtupleManager::CreateNtupleDColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<double>* column 
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
G4int G4CsvNtupleManager::CreateNtupleSColumn(G4int ntupleId,
                                              const G4String& name)
{
/*
  if ( vector ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple columns of vector are not supported." ;
    G4Exception("G4CsvAnalysisManager::CreateNtupleSColumn", 
                "Analysis_W002", JustWarning, description);
    return 0;
  }                
*/
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple S column", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleSColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvNtupleManager::CreateNtupleSColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<std::string>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<std::string>* column 
      = ntupleDescription->fNtuple->create_column<std::string>(name);  
    ntupleDescription->fNtupleSColumnMap[index] = column;
  }
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL2()->Message("create", "ntuple S column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
void G4CsvNtupleManager::FinishNtuple(G4int ntupleId)
{ 
  // nothing to be done if hippo header is inactivated
  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "FinishNtuple");
  if ( ! ntupleDescription )  return;

  // Write hippo header if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    if ( ! WriteHeader(ntupleDescription->fNtuple) ) {
       G4ExceptionDescription description;
       description << "      " 
                   << "Writing ntuple header has failed. ";
       G4Exception("G4CsvNtupleManager::Finish()",
                   "Analysis_W022", JustWarning, description);
    }
  }
}
     
//_____________________________________________________________________________
G4bool G4CsvNtupleManager::FillNtupleIColumn(G4int columnId, G4int value)
{
  return FillNtupleIColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::FillNtupleFColumn(G4int columnId, G4float value)
{
  return FillNtupleFColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::FillNtupleDColumn(G4int columnId, G4double value)
{
  return FillNtupleDColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::FillNtupleSColumn(G4int columnId, 
                                             const G4String& value)
{
  return FillNtupleSColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::AddNtupleRow()
{ 
  return AddNtupleRow(fFirstId);
}  

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::FillNtupleIColumn(G4int ntupleId, G4int columnId, 
                                               G4int value)
{
  tools::wcsv::ntuple::column<int>* column 
    = GetNtupleIColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::FillNtupleIColumn()",
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
G4bool G4CsvNtupleManager::FillNtupleFColumn(G4int ntupleId, G4int columnId, 
                                               G4float value)
{
  tools::wcsv::ntuple::column<float>* column 
    = GetNtupleFColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::FillNtupleFColumn()",
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
G4bool G4CsvNtupleManager::FillNtupleDColumn(G4int ntupleId, G4int columnId, 
                                               G4double value)
{
   tools::wcsv::ntuple::column<double>* column 
     = GetNtupleDColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::FillNtupleDColumn()",
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
G4bool G4CsvNtupleManager::FillNtupleSColumn(G4int ntupleId, G4int columnId, 
                                             const G4String& value)
{
   tools::wcsv::ntuple::column<std::string>* column 
     = GetNtupleSColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4CsvNtupleManager::FillNtupleSColumn()",
                "Analysis_W011", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fState.GetVerboseL4()->Message("fill", "ntuple S column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvNtupleManager::AddNtupleRow(G4int ntupleId)
{ 
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "AddNtupleRow");
  if ( ! ntupleDescription ) return false;

  if ( ! ntupleDescription->fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4CsvNtupleManager::AddNtupleRow()",
                "Analysis_W022", JustWarning, description);
    return false;
  }  
  
  ntupleDescription->fNtuple->add_row();
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description);
  }  
#endif

  return true;
}
