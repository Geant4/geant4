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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

G4ThreadLocal G4CsvAnalysisManager* G4CsvAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4CsvAnalysisManager* G4CsvAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new G4CsvAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager()
 : G4VAnalysisManager("Csv"),
   fNtupleVector()
{
  if ( fgInstance ) {
    G4ExceptionDescription description;
    description << "      " 
                << "G4CsvAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("G4CsvAnalysisManager::G4CsvAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  fgInstance = this;
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{  
  std::vector<G4CsvNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    delete (*it);
  }   

  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetNtupleFileName(
                                  G4CsvNtupleDescription* ntupleDescription) const
{                                  
  G4String name(fFileName);
  name.append("_");
  name.append(ntupleDescription->fNtupleBooking->m_name);
  // Add file extension .csv if no extension is given
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }
  return name;
}      

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CreateNtupleFile(
                                  G4CsvNtupleDescription* ntupleDescription)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "file", GetNtupleFileName(ntupleDescription));
#endif

  std::ofstream* ntupleFile 
    = new std::ofstream(GetNtupleFileName(ntupleDescription));
  if ( ntupleFile->fail() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " 
                << GetNtupleFileName(ntupleDescription);
    G4Exception("G4CsvAnalysisManager::CreateNtupleFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "file", GetNtupleFileName(ntupleDescription));
#endif

  ntupleDescription->fFile = ntupleFile;
  return true;
}  

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseNtupleFile(
                                  G4CsvNtupleDescription* ntupleDescription)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "file", GetNtupleFileName(ntupleDescription));
#endif

  // close file
  ntupleDescription->fFile->close(); 

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", GetNtupleFileName(ntupleDescription));
#endif

  return true; 
} 
   
//_____________________________________________________________________________
void G4CsvAnalysisManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( ! fNtupleVector.size() ) return;     
  
  std::vector<G4CsvNtupleDescription*>::iterator itn;  
  for (itn = fNtupleVector.begin(); itn != fNtupleVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;  
    if ( ! ntupleBooking ) continue;
    
#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif

    // create a file for this ntuple
    if ( ! CreateNtupleFile((*itn)) ) continue;

    // create ntuple
    (*itn)->fNtuple
      = new tools::wcsv::ntuple(*((*itn)->fFile), G4cerr, *ntupleBooking);
  
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
          G4Exception("G4CsvAnalysisManager::CreateNtupleFromBooking()",
                      "Analysis_W004", JustWarning, description);
        }
      }
    }
#ifdef G4VERBOSE
    if ( fpVerboseL3 ) 
      fpVerboseL3->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif
  }
}   

//_____________________________________________________________________________
tools::wcsv::ntuple::column<int>*    
G4CsvAnalysisManager::GetNtupleIColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
  if ( ! ntupleDecription ) {
    // add exception
    return 0;
  }    

  std::map<G4int, tools::wcsv::ntuple::column<int>* >& ntupleIColumnMap
    = ntupleDecription->fNtupleIColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<int>* >::const_iterator it
    = ntupleIColumnMap.find(columnId);
  if ( it == ntupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wcsv::ntuple::column<float>*  
G4CsvAnalysisManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleFColumn");
  if ( ! ntupleDecription ) {
    // add exception
    return 0;
  }    

  std::map<G4int, tools::wcsv::ntuple::column<float>* >& ntupleFColumnMap
    = ntupleDecription->fNtupleFColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<float>* >::const_iterator it
    = ntupleFColumnMap.find(columnId);
  if ( it == ntupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wcsv::ntuple::column<double>* 
G4CsvAnalysisManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  G4CsvNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleDColumn");
  if ( ! ntupleDecription ) {
    // add exception
    return 0;
  }    

  std::map<G4int, tools::wcsv::ntuple::column<double>* >& ntupleDColumnMap
    = ntupleDecription->fNtupleDColumnMap;
  std::map<G4int, tools::wcsv::ntuple::column<double>* >::const_iterator it
    = ntupleDColumnMap.find(columnId);
  if ( it == ntupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Reset()
{
  std::vector<G4CsvNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    delete (*it)->fNtuple; 
    (*it)->fNtuple = 0;
  }  
  
  return true;
}  
 
//_____________________________________________________________________________
void  G4CsvAnalysisManager::ExceptionForHistograms(
                                        const G4String& functionName) const
{
  G4String inFunction = "G4CsvAnalysisManager::";
  inFunction += functionName;
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception(inFunction, "Analysis_W005", JustWarning, description);
}  

//_____________________________________________________________________________
G4CsvNtupleDescription* G4CsvAnalysisManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool /*onlyIfActive*/) const
{                                      
  G4int index = id - fFirstNtupleId;
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4CsvAnalysisManager::";
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
G4bool G4CsvAnalysisManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// To be added: ntuple

  return true;
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;

  // Create ntuples if they are booked  
  if ( fNtupleVector.size() ) 
    CreateNtuplesFromBooking();
  
  fLockFileName = true;
  fLockNtupleDirectoryName = true;
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Write() 
{
  // nothing to be done for Csv file
  G4bool result = true;
  
  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii();
  }   

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFile()
{
  G4bool result = true;

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "files", "");
#endif

  // close files
  std::vector<G4CsvNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    CloseNtupleFile((*it));
  }  

  // reset data
  result = Reset();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Resetting data failed";
    G4Exception("G4CsvAnalysisManager::CloseFile()",
              "Analysis_W002", JustWarning, description);
    result = false;       
  } 

  fLockFileName = false;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    fpVerboseL1->Message("close", "files", "");
  }  
#endif

  return result; 
} 
   
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               const G4String& /*unitName*/, 
                               const G4String& /*fcnName*/)
{
  ExceptionForHistograms("CreateH1");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/)
{
  ExceptionForHistograms("CreateH2");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1(G4int /*id*/,
                                G4int /*nbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/,
                                const G4String& /*unitName*/, 
                                const G4String& /*fcnName*/)
{                                
  ExceptionForHistograms("SetH1");
  return false;
}
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2(G4int /*id*/,
                                G4int /*nxbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/, 
                                G4int /*nybins*/, 
                                G4double /*ymin*/, G4double /*ymax*/,
                                const G4String& /*xunitName*/, 
                                const G4String& /*yunitName*/, 
                                const G4String& /*xfcnName*/,
                                const G4String& /*yfcnName*/)
{                                
  ExceptionForHistograms("SetH2");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::ScaleH1(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH1");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::ScaleH2(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH2");
  return false;
}

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtuple(const G4String& name, 
                                         const G4String& title)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple", name);
#endif

  // Create ntuple description
  G4int index = fNtupleVector.size();
  G4CsvNtupleDescription* ntupleDescription
    = new G4CsvNtupleDescription();
  fNtupleVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking = new tools::ntuple_booking();
  ntupleDescription->fNtupleBooking->m_name = name;
  ntupleDescription->fNtupleBooking->m_title = title;
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open (what means here that
  // a filename was already set)
  if ( fFileName.size() ) {
    if ( CreateNtupleFile(ntupleDescription) ) {
      ntupleDescription->fNtuple 
        = new tools::wcsv::ntuple(*(ntupleDescription->fFile));
           // ntuple object is deleted when closing a file
    }       
  }  

  fLockFirstNtupleId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << index + fFirstNtupleId;
    fpVerboseL2->Message("create", "ntuple", description);
  } 
#endif

  return index + fFirstNtupleId;
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleIColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleFColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleDColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleDColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void G4CsvAnalysisManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleIColumn(G4int ntupleId, const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple I column", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleIColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleIColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<int>* column 
      = ntupleDescription->fNtuple->create_column<int>(name);  
    ntupleDescription->fNtupleIColumnMap[index] = column;
  }  

  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL2->Message("create", "ntuple I column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleFColumn(G4int ntupleId, const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple F column", description);
  } 
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleFColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleFColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<float>* column 
      = ntupleDescription->fNtuple->create_column<float>(name);  
    ntupleDescription->fNtupleFColumnMap[index] = column;
  }  

  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL2->Message("create", "ntuple F column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleDColumn(G4int ntupleId, const G4String& name)   
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple D column", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::wcsv::ntuple::column<double>* column 
      = ntupleDescription->fNtuple->create_column<double>(name);  
    ntupleDescription->fNtupleDColumnMap[index] = column;
  }
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL2->Message("create", "ntuple D column", description);
  }  
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
void G4CsvAnalysisManager::FinishNtuple(G4int /*ntupleId*/)
{ 
  // nothing to be done here
}
     
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillH1(G4int /*id*/, 
                                    G4double /*value*/, G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::FillH1()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillH2(G4int /*id*/, 
                                    G4double /*xvalue*/, G4double /*yvalue*/,
                                    G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::FillH2()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleIColumn(G4int columnId, G4int value)
{
  return FillNtupleIColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleFColumn(G4int columnId, G4float value)
{
  return FillNtupleFColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleDColumn(G4int columnId, G4double value)
{
  return FillNtupleDColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleIColumn(G4int ntupleId, G4int columnId, 
                                               G4int value)
{
  tools::wcsv::ntuple::column<int>* column 
    = GetNtupleIColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fpVerboseL4->Message("fill", "ntuple I column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleFColumn(G4int ntupleId, G4int columnId, 
                                               G4float value)
{
  tools::wcsv::ntuple::column<float>* column 
    = GetNtupleFColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fpVerboseL4->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleDColumn(G4int ntupleId, G4int columnId, 
                                               G4double value)
{
   tools::wcsv::ntuple::column<double>* column 
     = GetNtupleDColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "ntupleId " <<  ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId  
                << " columnId " << columnId << " value " << value;
    fpVerboseL4->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::AddNtupleRow()
{ 
  return AddNtupleRow(fFirstNtupleId);
}  

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::AddNtupleRow(G4int ntupleId)
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fpVerboseL4->Message("add", "ntuple row", description);
  }  
#endif

  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "AddNtupleRow");

  if ( ! ntupleDescription || ! ntupleDescription->fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4CsvAnalysisManager::AddNtupleRow()",
                "Analysis_W008", JustWarning, description);
    return false;
  }  
  
  ntupleDescription->fNtuple->add_row();
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fpVerboseL4->Message("add", "ntuple row", description);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
tools::wcsv::ntuple* G4CsvAnalysisManager::GetNtuple() const
{
  return GetNtuple(fFirstNtupleId);
}  

//_____________________________________________________________________________
tools::wcsv::ntuple* G4CsvAnalysisManager::GetNtuple(G4int ntupleId) const
{
  G4CsvNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");
    
  return ntupleDescription->fNtuple;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH1Nbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Nbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Width(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xwidth");
  return 0;
}
  
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH2Nxbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NXbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2XWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XWidth");
  return 0;
}
  
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH2Nybins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NYbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Ymin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Ymax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2YWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YWidth");
  return 0;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1Title(G4int /*id*/, 
                                        const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1XAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1YAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2Title(G4int /*id*/, 
                                        const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2XAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2YAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2ZAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2ZAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1YAxisTitle");
  return "";
}


//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2ZAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2ZAxisTitle");
  return "";
}
  
