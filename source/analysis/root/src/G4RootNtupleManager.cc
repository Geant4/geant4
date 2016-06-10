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
// $Id: G4RootNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootNtupleManager.hh"
#include "G4RootNtupleDescription.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state)
 : G4VNtupleManager(state),
   fNtupleDirectory(0),
   fNtupleDescriptionVector(),
   fNtupleVector()
{
}

//_____________________________________________________________________________
G4RootNtupleManager::~G4RootNtupleManager()
{  
  std::vector<G4RootNtupleDescription*>::iterator it;  
  for (it = fNtupleDescriptionVector.begin(); it != fNtupleDescriptionVector.end(); it++ ) {
    delete (*it);
  }   
}

// 
// private methods
//

//_____________________________________________________________________________
tools::wroot::ntuple::column<int>*    
G4RootNtupleManager::GetNtupleIColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wroot::ntuple::column<int>* >& ntupleIColumnMap
    = ntupleDecription->fNtupleIColumnMap;
  std::map<G4int, tools::wroot::ntuple::column<int>* >::const_iterator it
    = ntupleIColumnMap.find(columnId);

  if ( it == ntupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::GetNtupleIColumn()",
                "Analysis_W001", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wroot::ntuple::column<float>*  
G4RootNtupleManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleFColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wroot::ntuple::column<float>* >& ntupleFColumnMap
    = ntupleDecription->fNtupleFColumnMap;
  std::map<G4int, tools::wroot::ntuple::column<float>* >::const_iterator it
    = ntupleFColumnMap.find(columnId);
    
  if ( it == ntupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::GetNtupleFColumn()",
                "Analysis_W001", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wroot::ntuple::column<double>* 
G4RootNtupleManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleDColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wroot::ntuple::column<double>* >& ntupleDColumnMap
    = ntupleDecription->fNtupleDColumnMap;
  std::map<G4int, tools::wroot::ntuple::column<double>* >::const_iterator it
    = ntupleDColumnMap.find(columnId);
    
  if ( it == ntupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::GetNtupleDColumn()",
                "Analysis_W001", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
tools::wroot::ntuple::column_string* 
G4RootNtupleManager::GetNtupleSColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleSColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::wroot::ntuple::column_string* >& ntupleSColumnMap
    = ntupleDecription->fNtupleSColumnMap;
  std::map<G4int, tools::wroot::ntuple::column_string* >::const_iterator it
    = ntupleSColumnMap.find(columnId);
    
  if ( it == ntupleSColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::GetNtupleSColumn()",
                "Analysis_W001", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
G4RootNtupleDescription* G4RootNtupleManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool /*onlyIfActive*/) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootNtupleManager::";
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
// protected methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( ! fNtupleDescriptionVector.size() ) return;     
  
  std::vector<G4RootNtupleDescription*>::iterator itn;  
  for (itn = fNtupleDescriptionVector.begin(); itn != fNtupleDescriptionVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;
    if ( ! ntupleBooking ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "ntuple", ntupleBooking->name());
#endif

    (*itn)->fNtuple
      = new tools::wroot::ntuple(*fNtupleDirectory, *ntupleBooking);
    fNtupleVector.push_back((*itn)->fNtuple);  
  
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
        else if ( it->cls_id()== tools::_cid(double(0))) {
          (*itn)->fNtupleDColumnMap[index++] 
            = (*itn)->fNtuple->find_column<double>(it->name());
        }
        else if ( it->cls_id()== tools::_cid(std::string(""))) {
          (*itn)->fNtupleSColumnMap[index++] 
            = (*itn)->fNtuple->find_column_string(it->name());
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << it->name();
          G4Exception("G4RootNtupleManager::CreateNtupleFromBooking()",
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
G4bool G4RootNtupleManager::IsEmpty() const
{
  return ! fNtupleDescriptionVector.size();
}  
 
//_____________________________________________________________________________
G4bool G4RootNtupleManager::Reset()
{
// Reset ntuples

  std::vector<G4RootNtupleDescription*>::iterator it;  
  for (it = fNtupleDescriptionVector.begin(); it != fNtupleDescriptionVector.end(); it++ ) {
    // ntuple is deleted automatically when file is closed
    // delete (*it)->fNtuple;
    (*it)->fNtuple=0; 
  } 
  fNtupleVector.clear(); 
  
  return true;
}  
 
//_____________________________________________________________________________
tools::wroot::ntuple* G4RootNtupleManager::GetNtuple() const
{
  return GetNtuple(fFirstId);
}  

//_____________________________________________________________________________
tools::wroot::ntuple* G4RootNtupleManager::GetNtuple(G4int ntupleId) const
{
  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");
    
  if ( ! ntupleDescription ) return 0; 
    
  return ntupleDescription->fNtuple;  
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtuple(const G4String& name, 
                                        const G4String& title)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple", name);
#endif

  // Create ntuple description
  G4int index = fNtupleDescriptionVector.size();
  G4RootNtupleDescription* ntupleDescription
    = new G4RootNtupleDescription();
  fNtupleDescriptionVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking 
    = new tools::ntuple_booking(name, title);
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open
  if ( fNtupleDirectory ) {
    ntupleDescription->fNtuple 
      = new tools::wroot::ntuple(*fNtupleDirectory, name, title);
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);       
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
G4int G4RootNtupleManager::CreateNtupleIColumn(const G4String& name,
                                               std::vector<int>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleIColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleFColumn(const G4String& name,
                                               std::vector<float>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleFColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleDColumn(const G4String& name,
                                               std::vector<double>* vector)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleDColumn(ntupleId, name, vector);
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleSColumn(const G4String& name)
{
  G4int ntupleId = fNtupleDescriptionVector.size() + fFirstId - 1;
  return CreateNtupleSColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void G4RootNtupleManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleIColumn(G4int ntupleId, 
                                               const G4String& name,
                                               std::vector<int>* vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple I column", description);
  }  
#endif
  
  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleIColumn");
  if ( ! ntupleDescription )  return kInvalidId;   
    
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;
  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleIColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  if ( ! vector )
    ntupleBooking->add_column<int>(name);
  else
    ntupleBooking->add_column<int>(name, *vector);
      
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    if ( ! vector ) {
      tools::wroot::ntuple::column<int>* column
        = ntupleDescription->fNtuple->create_column<int>(name);
      ntupleDescription->fNtupleIColumnMap[index] = column;
    }  
    else {
      ntupleDescription->fNtuple->create_column<int>(name, *vector); 
    }   
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
G4int G4RootNtupleManager::CreateNtupleFColumn(G4int ntupleId, 
                                               const G4String& name,
                                               std::vector<float>* vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple F column", description);
  } 
#endif

  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleFColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleFColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  if ( ! vector )
    ntupleBooking->add_column<float>(name);
  else
    ntupleBooking->add_column<float>(name, *vector);
      
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    if ( ! vector ) {
      tools::wroot::ntuple::column<float>* column
        = ntupleDescription->fNtuple->create_column<float>(name);
      ntupleDescription->fNtupleFColumnMap[index] = column;
    }  
    else {
      ntupleDescription->fNtuple->create_column<float>(name, *vector); 
    }   
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
G4int G4RootNtupleManager::CreateNtupleDColumn(G4int ntupleId, 
                                               const G4String& name,  
                                               std::vector<double>* vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple D column", description);
  }  
#endif

  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleDColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  if ( ! vector )
    ntupleBooking->add_column<double>(name);  
  else
    ntupleBooking->add_column<double>(name, *vector);  
    
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    if ( ! vector ) {
      tools::wroot::ntuple::column<double>* column
        = ntupleDescription->fNtuple->create_column<double>(name);
      ntupleDescription->fNtupleDColumnMap[index] = column;
    }  
    else {
      ntupleDescription->fNtuple->create_column<double>(name, *vector); 
    }   
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
G4int G4RootNtupleManager::CreateNtupleSColumn(G4int ntupleId, 
                                               const G4String& name)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fState.GetVerboseL4()->Message("create", "ntuple S column", description);
  }  
#endif

  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleSColumn");
  if ( ! ntupleDescription )  return kInvalidId;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleSColumn()",
                "Analysis_W002", JustWarning, description);
    return kInvalidId;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->columns().size();
  ntupleBooking->add_column<std::string>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wroot::ntuple::column_string* column
      = ntupleDescription->fNtuple->create_column_string(name);
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
void G4RootNtupleManager::FinishNtuple(G4int /*ntupleId*/)
{ 
  // nothing to be done here
}
   
//_____________________________________________________________________________
G4bool G4RootNtupleManager::FillNtupleIColumn(G4int columnId, G4int value)
{
  return FillNtupleIColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootNtupleManager::FillNtupleFColumn(G4int columnId, G4float value)
{
  return FillNtupleFColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootNtupleManager::FillNtupleDColumn(G4int columnId, G4double value)
{
  return FillNtupleDColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootNtupleManager::FillNtupleSColumn(G4int columnId, 
                                              const G4String& value)
{
  return FillNtupleSColumn(fFirstId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootNtupleManager::AddNtupleRow()
{ 
  return AddNtupleRow(fFirstId);
}  

//_____________________________________________________________________________
G4bool G4RootNtupleManager::FillNtupleIColumn(G4int ntupleId, G4int columnId, 
                                                G4int value)
{
  tools::wroot::ntuple::column<int>* column 
    = GetNtupleIColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::FillNtupleIColumn()",
                "Analysis_W001", JustWarning, description);
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
G4bool G4RootNtupleManager::FillNtupleFColumn(G4int ntupleId, G4int columnId, 
                                                G4float value)
{
  tools::wroot::ntuple::column<float>* column 
    = GetNtupleFColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::FillNtupleFColumn()",
                "Analysis_W001", JustWarning, description);
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
G4bool G4RootNtupleManager::FillNtupleDColumn(G4int ntupleId, G4int columnId, 
                                                G4double value)
{
  tools::wroot::ntuple::column<double>* column 
    = GetNtupleDColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::FillNtupleDColumn()",
                "Analysis_W001", JustWarning, description);
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
G4bool G4RootNtupleManager::FillNtupleSColumn(G4int ntupleId, G4int columnId, 
                                                const G4String& value)
{
  tools::wroot::ntuple::column_string* column 
    = GetNtupleSColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << " columnId " << columnId << " does not exist.";
    G4Exception("G4RootNtupleManager::FillNtupleSColumn()",
                "Analysis_W001", JustWarning, description);
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
G4bool G4RootNtupleManager::AddNtupleRow(G4int ntupleId)
{ 
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description);
  }  
#endif

  G4RootNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "AddNtupleRow");
  if ( ! ntupleDescription ) return false;

  if ( ! ntupleDescription->fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << " does not exist. ";
    G4Exception("G4RootNtupleManager::AddNtupleRow()",
                "Analysis_W022", JustWarning, description);
    return false;
  }  
  
  G4bool result = ntupleDescription->fNtuple->add_row();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << "adding row has failed.";
    G4Exception("G4RootNtupleManager::AddNtupleRow()",
                "Analysis_W002", JustWarning, description);
  }         
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "ntuple row", description, result);
  }  
#endif

  return result;
}
