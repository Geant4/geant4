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

#include "tools/wroot/file"

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state)
 : G4VNtupleManager(state),
   fNtupleDirectory(0),
   fNtupleVector()
{
}

//_____________________________________________________________________________
G4RootNtupleManager::~G4RootNtupleManager()
{  
  std::vector<G4RootNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
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
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wroot::ntuple::column<float>*  
G4RootNtupleManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
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
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wroot::ntuple::column<double>* 
G4RootNtupleManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  G4RootNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
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
                "Analysis_W009", JustWarning, description);
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
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  return fNtupleVector[index];
}
  
// 
// protected methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( ! fNtupleVector.size() ) return;     
  
  std::vector<G4RootNtupleDescription*>::iterator itn;  
  for (itn = fNtupleVector.begin(); itn != fNtupleVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;
    if ( ! ntupleBooking ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif

    (*itn)->fNtuple
      = new tools::wroot::ntuple(*fNtupleDirectory, *ntupleBooking);
  
    if ( ntupleBooking->m_columns.size() ) {
      // store ntuple columns in local maps
      const std::vector<tools::ntuple_booking::col_t>& columns 
        = ntupleBooking->m_columns;
      std::vector<tools::ntuple_booking::col_t>::const_iterator it;
      G4int index = 0;
      for ( it = columns.begin(); it!=columns.end(); ++it) {
        if ( (*it).second == tools::_cid(int(0) ) ) {
          (*itn)->fNtupleIColumnMap[index++] 
            = (*itn)->fNtuple->find_column<int>((*it).first);
        }
        else if ( (*it).second == tools::_cid(float(0) ) ) {
          (*itn)->fNtupleFColumnMap[index++] 
            = (*itn)->fNtuple->find_column<float>((*it).first);
        } 
        else if ( (*it).second== tools::_cid(double(0))) {
          (*itn)->fNtupleDColumnMap[index++] 
            = (*itn)->fNtuple->find_column<double>((*it).first);
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << (*it).first;
          G4Exception("G4RootNtupleManager::CreateNtupleFromBooking()",
                      "Analysis_W004", JustWarning, description);
        }
      }
    }
#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()
        ->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif
  }
}   

//_____________________________________________________________________________
G4bool G4RootNtupleManager::IsEmpty() const
{
  return ! fNtupleVector.size();
}  
 
//_____________________________________________________________________________
G4bool G4RootNtupleManager::Reset()
{
// Reset ntuples

  std::vector<G4RootNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    // ntuple is deleted automatically when file is closed
    // delete (*it)->fNtuple;
    (*it)->fNtuple=0; 
  }  
  
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
  G4int index = fNtupleVector.size();
  G4RootNtupleDescription* ntupleDescription
    = new G4RootNtupleDescription();
  fNtupleVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking = new tools::ntuple_booking();
  ntupleDescription->fNtupleBooking->m_name = name;
  ntupleDescription->fNtupleBooking->m_title = title;
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open
  if ( fNtupleDirectory ) {
    ntupleDescription->fNtuple 
      = new tools::wroot::ntuple(*fNtupleDirectory, name, title);
           // ntuple object is deleted automatically when closing a file
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
G4int G4RootNtupleManager::CreateNtupleIColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstId - 1;
  return CreateNtupleIColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleFColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstId - 1;
  return CreateNtupleFColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleDColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstId - 1;
  return CreateNtupleDColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void G4RootNtupleManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
//_____________________________________________________________________________
G4int G4RootNtupleManager::CreateNtupleIColumn(G4int ntupleId, 
                                                 const G4String& name)
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
  if ( ! ntupleDescription )  return -1;   
    
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;
  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleIColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wroot::ntuple::column<int>* column 
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
G4int G4RootNtupleManager::CreateNtupleFColumn(G4int ntupleId, const G4String& name)
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
  if ( ! ntupleDescription )  return -1;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleFColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wroot::ntuple::column<float>* column 
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
G4int G4RootNtupleManager::CreateNtupleDColumn(G4int ntupleId, const G4String& name)   
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
  if ( ! ntupleDescription )  return -1;   

  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4RootNtupleManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wroot::ntuple::column<double>* column 
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
                "Analysis_W009", JustWarning, description);
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
                "Analysis_W009", JustWarning, description);
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
                "Analysis_W009", JustWarning, description);
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
                "Analysis_W008", JustWarning, description);
    return false;
  }  
  
  G4bool result = ntupleDescription->fNtuple->add_row();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << "adding row has failed.";
    G4Exception("G4RootNtupleManager::AddNtupleRow()",
                "Analysis_W004", JustWarning, description);
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
