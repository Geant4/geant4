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

#include "G4XmlAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "tools/waxml/begend"
#include "tools/waxml/histos"

#include <iostream>

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH2Mutex = G4MUTEX_INITIALIZER;
}  

G4XmlAnalysisManager* G4XmlAnalysisManager::fgMasterInstance = 0;
G4ThreadLocal G4XmlAnalysisManager* G4XmlAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4XmlAnalysisManager* G4XmlAnalysisManager::Create(G4bool isMaster)
{
  if ( fgInstance == 0 ) {
    fgInstance = new G4XmlAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4XmlAnalysisManager* G4XmlAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new G4XmlAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4XmlAnalysisManager::G4XmlAnalysisManager(G4bool isMaster)
 : G4VAnalysisManager("Xml"),
   fIsMaster(isMaster),
   fHnFile(0),
   fH1Vector(),   
   fH2Vector(),   
   fH1NameIdMap(),  
   fH2NameIdMap(),  
   fNtupleVector()
{
  if ( ( isMaster && fgMasterInstance ) ||
       ( (! isMaster ) && fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4XmlAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4XmlAnalysisManager::G4XmlAnalysisManager",
                "Analysis_F001", FatalException, description);
  }              
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;
}

//_____________________________________________________________________________
G4XmlAnalysisManager::~G4XmlAnalysisManager()
{  
  std::vector<tools::histo::h1d*>::iterator it;
  for ( it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  } 
   
  std::vector<tools::histo::h2d*>::iterator it2;
  for ( it2 = fH2Vector.begin(); it2 != fH2Vector.end(); it2++ ) {
    delete *it2;
  }
  
  std::vector<G4XmlNtupleDescription*>::iterator it3;  
  for (it3 = fNtupleVector.begin(); it3 != fNtupleVector.end(); it3++ ) {
    delete (*it3);
  }   
   
  delete fHnFile;  

  if ( fIsMaster ) fgMasterInstance = 0;
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetHnFileName() const
{                                  
  G4String name(fFileName);
  // Add thread Id to a file name if MT processing
  if ( ! fIsMaster ) {
    std::ostringstream os;
    os << G4GetPidId();
    name.append("_t");
    name.append(os.str());
  }  
  // Add file extension .xml if no extension is given
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }
  return name;
}      

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CreateHnFile()
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "file", GetHnFileName());
#endif
  
  // delete a previous file if it exists
  if ( fHnFile ) delete fHnFile; 
  
  fHnFile = new std::ofstream(GetHnFileName());
  if ( fHnFile->fail() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << GetHnFileName();
    G4Exception("G4XmlAnalysisManager::OpenFile()",
              "Analysis_W001", JustWarning, description);
    return false;
  }

  tools::waxml::begin(*fHnFile);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "file", GetHnFileName());
#endif

  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseHnFile()
{
  // No file may be open if no master manager is instantiated
  // and no histograms were booked
  if ( ! fHnFile ) return true;

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "file", GetHnFileName());
#endif

  // close file
  tools::waxml::end(*fHnFile);
  fHnFile->close(); 

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", GetHnFileName());
#endif

  return true; 
} 
   
//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetNtupleFileName(
                                    G4XmlNtupleDescription* ntupleDescription) const
{                                  
  G4String name(fFileName);
  name.append("_");
  name.append(ntupleDescription->fNtupleBooking->m_name);
  // Add thread Id to a file name if MT processing
  if ( ! fIsMaster ) {
    std::ostringstream os;
    os << G4GetPidId();
    name.append("_t");
    name.append(os.str());
  }  
  // Add file extension .xml if no extension is given
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }
  return name;
}      

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CreateNtupleFile(
                                    G4XmlNtupleDescription* ntupleDescription)
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
    G4Exception("G4XmlAnalysisManager::CreateNtupleFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  
  tools::waxml::begin(*ntupleFile);
  ntupleDescription->fFile = ntupleFile;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "file", GetNtupleFileName(ntupleDescription));
#endif

  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseNtupleFile(
                                    G4XmlNtupleDescription* ntupleDescription)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "file", GetNtupleFileName(ntupleDescription));
#endif

  // close file
  tools::waxml::end(*(ntupleDescription->fFile));
  ntupleDescription->fFile->close(); 

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", GetNtupleFileName(ntupleDescription));
#endif

  return true; 
} 
   
//_____________________________________________________________________________
void G4XmlAnalysisManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( ! fNtupleVector.size() ) return;     
  
  std::vector<G4XmlNtupleDescription*>::iterator itn;  
  for (itn = fNtupleVector.begin(); itn != fNtupleVector.end(); itn++ ) {

    tools::ntuple_booking* ntupleBooking = (*itn)->fNtupleBooking;
    if ( ! ntupleBooking ) continue;

// Create ntuple from ntuple_booking.
#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif

    // create a file for this ntuple
    if ( ! CreateNtupleFile((*itn)) ) continue;

    // create ntuple
    (*itn)->fNtuple
      = new tools::waxml::ntuple(*((*itn)->fFile), G4cerr, *ntupleBooking);

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
        else if( (*it).second == tools::_cid(float(0) ) ) {
          (*itn)->fNtupleFColumnMap[index++] 
            = (*itn)->fNtuple->find_column<float>((*it).first);
        } 
        else if((*it).second== tools::_cid(double(0))) {
          (*itn)->fNtupleDColumnMap[index++] 
            = (*itn)->fNtuple->find_column<double>((*it).first);
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << (*it).first;
          G4Exception("G4XmlAnalysisManager::CreateNtuplesFromBooking()",
                    "Analysis_W004", JustWarning, description);
        }
      }
    }
    FinishNtuple();
#ifdef G4VERBOSE
    if ( fpVerboseL3 ) 
      fpVerboseL3->Message("create from booking", "ntuple", ntupleBooking->m_name);
#endif
  }  
}   

//_____________________________________________________________________________
void G4XmlAnalysisManager::AddH1Vector(
                               std::vector<tools::histo::h1d*>& h1Vector)
{
#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("merge", "all h1", "");
#endif
  std::vector<tools::histo::h1d*>::iterator itw = h1Vector.begin();
  std::vector<tools::histo::h1d*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fpVerboseL1 ) 
      fpVerboseL1->Message("merge", "all h1", "");
#endif
}  

//_____________________________________________________________________________
void G4XmlAnalysisManager::AddH2Vector(
                               std::vector<tools::histo::h2d*>& h2Vector)
{
#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("merge", "all h2", "");
#endif
  std::vector<tools::histo::h2d*>::iterator itw = h2Vector.begin();
  std::vector<tools::histo::h2d*>::iterator it;
  for (it = fH2Vector.begin(); it != fH2Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fpVerboseL1 ) 
      fpVerboseL1->Message("merge", "all h2", "");
#endif
}  

//_____________________________________________________________________________
tools::waxml::ntuple::column<int>*    
G4XmlAnalysisManager::GetNtupleIColumn(G4int ntupleId, G4int columnId) const
{
  G4XmlNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleIColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::waxml::ntuple::column<int>* >& ntupleIColumnMap
    = ntupleDecription->fNtupleIColumnMap;
  std::map<G4int, tools::waxml::ntuple::column<int>* >::const_iterator it
    = ntupleIColumnMap.find(columnId);
  if ( it == ntupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::GetNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::waxml::ntuple::column<float>*  
G4XmlAnalysisManager::GetNtupleFColumn(G4int ntupleId, G4int columnId) const
{
  G4XmlNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleFColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::waxml::ntuple::column<float>* >& ntupleFColumnMap
    = ntupleDecription->fNtupleFColumnMap;
  std::map<G4int, tools::waxml::ntuple::column<float>* >::const_iterator it
    = ntupleFColumnMap.find(columnId);
  if ( it == ntupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  

//_____________________________________________________________________________
tools::waxml::ntuple::column<double>* 
G4XmlAnalysisManager::GetNtupleDColumn(G4int ntupleId, G4int columnId) const
{
  G4XmlNtupleDescription* ntupleDecription
    = GetNtupleInFunction(ntupleId, "GetNtupleDColumn");
  if ( ! ntupleDecription ) return 0;

  std::map<G4int, tools::waxml::ntuple::column<double>* >& ntupleDColumnMap
    = ntupleDecription->fNtupleDColumnMap;
  std::map<G4int, tools::waxml::ntuple::column<double>* >::const_iterator it
    = ntupleDColumnMap.find(columnId);
  if ( it == ntupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      "  << "ntupleId " << ntupleId
                << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::h1d*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  
  
  std::vector<tools::histo::h2d*>::iterator it2;
  for (it2 = fH2Vector.begin(); it2 != fH2Vector.end(); it2++ ) {
    G4bool result = (*it2)->reset();
    if ( ! result ) finalResult = false;
  }  

  std::vector<G4XmlNtupleDescription*>::iterator it3;  
  for (it3 = fNtupleVector.begin(); it3 != fNtupleVector.end(); it3++ ) {
    delete (*it3)->fNtuple;
    (*it3)->fNtuple = 0;
  }  
  
  return finalResult;
}  
 
//_____________________________________________________________________________
tools::histo::h1d*  G4XmlAnalysisManager::GetH1InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4XmlAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  // Do not return histogram if inactive 
  if ( fActivation && onlyIfActive && ( ! GetActivation(kH1, id) ) ) {
    return 0; 
  }  
  
  return fH1Vector[index];
}  
                                      
//_____________________________________________________________________________
tools::histo::h2d*  G4XmlAnalysisManager::GetH2InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4XmlAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }

  // Do not return histogram if inactive 
  if ( fActivation  && onlyIfActive && ( ! GetActivation(kH2, id) ) ) {
    return 0; 
  }  
  
  return fH2Vector[index];
}
  
//_____________________________________________________________________________
G4XmlNtupleDescription* G4XmlAnalysisManager::GetNtupleInFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool /*onlyIfActive*/) const
{                                      
  G4int index = id - fFirstNtupleId;
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4XmlAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  return fNtupleVector[index];
}
  
//_____________________________________________________________________________
void G4XmlAnalysisManager::UpdateTitle(G4String& title, 
                                        const G4String& unitName, 
                                        const G4String& fcnName) const
{
  if ( fcnName != "none" )  { title += " "; title += fcnName; title += "("; }
  if ( unitName != "none" ) { title += " ["; title += unitName; title += "]";}
  if ( fcnName != "none" )  { title += ")"; }
}  
                                                          
//
// protected methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// (Only H1 implemented by now)
// According to the implementation by Michel Maire, originally in
// extended examples.

  // h1 histograms
  for ( G4int i=0; i<G4int(fH1Vector.size()); ++i ) {
    G4int id = i + fFirstHistoId;
    G4HnInformation* info = GetH1Information(id); 
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->fAscii ) continue; 
    tools::histo::h1d* h1 = fH1Vector[i];

#ifdef G4VERBOSE
    if ( fpVerboseL3 ) 
      fpVerboseL3->Message("write on ascii", "h1d", info->fName);
#endif
  
    output << "\n  1D histogram " << id << ": " << h1->title() 
           << "\n \n \t     X \t\t     Y" << G4endl;
    
    for (G4int j=0; j< G4int(h1->axis().bins()); ++j) {
       output << "  " << j << "\t" 
              << h1->axis().bin_center(j) << "\t"
              << h1->bin_height(j) << G4endl;
    } 
  }
  
  return true;
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;

  // Create file for histograms
  G4bool result = true;
  if ( fIsMaster ) 
    result = CreateHnFile();

  // Create ntuples if they are booked
  if ( fNtupleVector.size() ) 
    CreateNtuplesFromBooking();

  fLockFileName = true;

  return result;
}  
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::Write() 
{
  // ntuples 
  for ( G4int i=0; i<G4int(fNtupleVector.size()); ++i ) {
    if ( fNtupleVector[i]->fNtuple ) fNtupleVector[i]->fNtuple->write_trailer();
  }

  if ( fIsMaster || ( ! fgMasterInstance ) )  {

    if ( ! fgMasterInstance && 
         ( fH1Vector.size() || fH2Vector.size() ) ) {

      G4ExceptionDescription description;
      description 
        << "      " << "No master G4RootAnalysisManager instance exists." 
        << G4endl 
        << "      " << "Histogram data will not be merged.";
        G4Exception("G4RootAnalysisManager::Write()",
                  "Analysis_W014", JustWarning, description);
                  
      // Create Hn file per thread
      G4bool result = CreateHnFile();
      if ( ! result ) return false;       
    }

    // h1 histograms
    for ( G4int i=0; i<G4int(fH1Vector.size()); ++i ) {
      G4int id = i + fFirstHistoId;
      G4HnInformation* info = GetH1Information(id); 
      // skip writing if activation is enabled and H1 is inactivated
      if ( fActivation && ( ! info->fActivation ) ) continue; 
      tools::histo::h1d* h1 = fH1Vector[i];
#ifdef G4VERBOSE
      if ( fpVerboseL3 ) 
        fpVerboseL3->Message("write", "h1d", info->fName);
#endif
      G4String path = "/";
      path.append(fHistoDirectoryName);
      G4bool result
        = tools::waxml::write(*fHnFile, *h1, path, info->fName);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << info->fName << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W003", JustWarning, description);
        return false;       
      } 
      fLockHistoDirectoryName = true;
    }
 
    // h2 histograms
    for ( G4int i=0; i<G4int(fH2Vector.size()); ++i ) {
      G4int id = i + fFirstHistoId;
      G4HnInformation* info = GetH2Information(id); 
      // skip writing if inactivated
      if ( fActivation && ( ! info->fActivation ) ) continue;
      tools::histo::h2d* h2 = fH2Vector[i];
 #ifdef G4VERBOSE
      if ( fpVerboseL3 ) 
        fpVerboseL3->Message("write", "h2d", info->fName);
#endif
      G4String path = "/";
      path.append(fHistoDirectoryName);
      G4bool result
        = tools::waxml::write(*fHnFile, *h2, path, info->fName);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << info->fName << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W003", JustWarning, description);
        return false;       
      } 
      fLockHistoDirectoryName = true;
    }
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH1(&mergeH1Mutex);
    fgMasterInstance->AddH1Vector(fH1Vector);
    lH1.unlock();
    
    G4AutoLock lH2(&mergeH2Mutex);
    fgMasterInstance->AddH2Vector(fH2Vector);
    lH2.unlock();
  }  
  G4bool result = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("write", "file", GetFullFileName(), result);
#endif

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii();
  }   

  return result;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseFile()
{
  G4bool result = true;

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "files", "");
#endif

  // close histogram file
  CloseHnFile();
  
  // close ntuple foles
  std::vector<G4XmlNtupleDescription*>::iterator it;  
  for (it = fNtupleVector.begin(); it != fNtupleVector.end(); it++ ) {
    CloseNtupleFile((*it));
  }  

  // reset data
  result = Reset();
  if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4XmlAnalysisManager::CloseFile()",
                "Analysis_W002", JustWarning, description);
      result = false;       
  } 

  fLockFileName = false;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "files", "");
#endif

  return result; 
} 
   
//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateH1(const G4String& name, const G4String& title, 
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "H1", name);
#endif
  G4int index = fH1Vector.size();
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  tools::histo::h1d* h1 
    = new tools::histo::h1d(title, nbins, fcn(xmin), fcn(xmax));
            // h1 objects are deleted in destructor and reset when 
            // closing a file.

  G4String axisTitle;
  UpdateTitle(axisTitle,unitName, fcnName);        
  h1->add_annotation(tools::histo::key_axis_x_title(), axisTitle);
             
  fH1Vector.push_back(h1);
  AddH1Information(name, unitName, fcnName, unit, fcn);

  fLockFirstHistoId = true;
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "H1", name);
#endif
  fH1NameIdMap[name] = index + fFirstHistoId;
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateH2(const G4String& name, const G4String& title, 
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "H2", name);
#endif
  G4int index = fH2Vector.size();
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  tools::histo::h2d* h2 
    = new tools::histo::h2d(title, 
                            nxbins, xfcn(xmin), xfcn(xmax), 
                            nybins, yfcn(ymin), yfcn(ymax));
            // h1 objects are deleted in destructor and reset when 
            // closing a file.

  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  h2->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
             
  fH2Vector.push_back(h2);
  AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                   xunit, yunit, xfcn, yfcn);

  fLockFirstHistoId = true;
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "H2", name);
#endif
  fH2NameIdMap[name] = index + fFirstHistoId;
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH1(G4int id,
                                G4int nbins, G4double xmin, G4double xmax,
                                const G4String& unitName, const G4String& fcnName)
{                                

  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1", false, false);
  if ( ! h1d ) return false;

  G4HnInformation* info = GetH1Information(id);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("configure", "H1", info->fName);
#endif

  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  h1d->configure(nbins, fcn(xmin), fcn(xmax));
  info->fXUnitName = unitName;
  info->fYUnitName = unitName;
  info->fXFcnName = fcnName;
  info->fYFcnName = fcnName;
  info->fXUnit = unit;
  info->fYUnit = unit;
  info->fXFcn = fcn;
  info->fYFcn = fcn;
  SetActivation(kH1, id, true); 

  G4String axisTitle;
  UpdateTitle(axisTitle,unitName, fcnName);        
  h1d->add_annotation(tools::histo::key_axis_x_title(), axisTitle);

  return true;
}
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH2(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{                                
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  G4HnInformation* info = GetH2Information(id);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("configure", "H2", info->fName);
#endif

  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  h2d->configure(nxbins, xfcn(xmin), xfcn(xmax), 
                 nybins, yfcn(ymin), yfcn(ymax));
  info->fXUnitName = xunitName;
  info->fYUnitName = yunitName;
  info->fXFcnName = xfcnName;
  info->fYFcnName = yfcnName;
  info->fXUnit = xunit;
  info->fYUnit = yunit;
  info->fXFcn = xfcn;
  info->fYFcn = yfcn;
  SetActivation(kH2, id, true); 
  
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  h2d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::ScaleH1(G4int id, G4double factor)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "ScaleH1", false, false);
  if ( ! h1d ) return false;

  return h1d->scale(factor);
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::ScaleH2(G4int id, G4double factor)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "ScaleH2", false, false);
  if ( ! h2d ) return false;
  
  return h2d->scale(factor);
}  
                           
//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtuple(const G4String& name, 
                                         const G4String& title)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple", name);
#endif

  // Create ntuple description
  G4int index = fNtupleVector.size();
  G4XmlNtupleDescription* ntupleDescription
    = new G4XmlNtupleDescription();
  fNtupleVector.push_back(ntupleDescription);  

  // Create ntuple booking
  ntupleDescription->fNtupleBooking = new tools::ntuple_booking();
  ntupleDescription->fNtupleBooking->m_name = name;
  ntupleDescription->fNtupleBooking->m_title = title;

  // Create ntuple if the file is open (what means here that
  // a filename was already set)
  if ( fFileName.size() ) {
    if ( CreateNtupleFile(ntupleDescription) ) {
      ntupleDescription->fNtuple 
        = new tools::waxml::ntuple(*(ntupleDescription->fFile));
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
G4int G4XmlAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleIColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleFColumn(ntupleId, name);
}  

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleDColumn(const G4String& name)
{
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  return CreateNtupleDColumn(ntupleId, name);
}  

//_____________________________________________________________________________
void G4XmlAnalysisManager::FinishNtuple()
{ 
  G4int ntupleId = fNtupleVector.size() + fFirstNtupleId - 1;
  FinishNtuple(ntupleId);
}
   
//_____________________________________________________________________________
G4int G4XmlAnalysisManager::CreateNtupleIColumn(G4int ntupleId, const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple I column", description);
  }  
 #endif

  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleIColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4XmlAnalysisManager::CreateNtupleIColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::waxml::ntuple::column<int>* column 
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
G4int G4XmlAnalysisManager::CreateNtupleFColumn(G4int ntupleId, const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple F column", description);
  } 
#endif

  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleFColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4XmlAnalysisManager::CreateNtupleFColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::waxml::ntuple::column<float>* column 
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
G4int G4XmlAnalysisManager::CreateNtupleDColumn(G4int ntupleId, const G4String& name)   
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("create", "ntuple D column", description);
  }  
#endif

  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "CreateNtupleDColumn");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      "
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4XmlAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = ntupleBooking->m_columns.size();
  ntupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( ntupleDescription->fNtuple ) {
    tools::waxml::ntuple::column<double>* column 
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
void G4XmlAnalysisManager::FinishNtuple(G4int ntupleId)
{ 
  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "FinishNtuple");
  tools::ntuple_booking* ntupleBooking
    = ntupleDescription->fNtupleBooking;  

  if ( ! ntupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple " << ntupleId << " has to be created first. ";
    G4Exception("G4XmlAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return;       
  }

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << ntupleBooking->m_name << " ntupleId " << ntupleId; 
    fpVerboseL4->Message("finish", "ntuple", description);
  }  
#endif

  G4String path = "/";
  path.append(fNtupleDirectoryName);
  ntupleDescription->fNtuple
    ->write_header(path, ntupleBooking->m_name, ntupleBooking->m_title);  

  fLockNtupleDirectoryName = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << ntupleBooking->m_name << " ntupleId " << ntupleId; 
    fpVerboseL2->Message("finish", "ntuple", description);
  }  
#endif
}
   
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillH1(G4int id, G4double value, G4double weight)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "FillH1", true, false);
  if ( ! h1d ) return false;

  if ( fActivation && ( ! GetActivation(kH1, id) ) ) {
    //G4cout << "Skipping FillH1 for " << id << G4endl; 
    return false; 
  }  

  G4HnInformation* info = GetInformation(kH1, id);
  h1d->fill(info->fXFcn(value/info->fXUnit), weight);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value/GetXUnit(kH1, id);
    fpVerboseL4->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillH2(G4int id, 
                                    G4double xvalue, G4double yvalue, 
                                    G4double weight)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "FillH2", true, false);
  if ( ! h2d ) return false;

  if ( fActivation && ( ! GetActivation(kH2, id) ) ) return false; 

  G4HnInformation* info = GetInformation(kH2, id);
  h2d->fill(info->fXFcn(xvalue/info->fXUnit), 
            info->fYFcn(yvalue/info->fYUnit), weight);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue/GetXUnit(kH2, id) 
                << " yvalue " << yvalue/GetYUnit(kH2, id);
    fpVerboseL4->Message("fill", "H2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleIColumn(G4int columnId, G4int value)
{
  return FillNtupleIColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleFColumn(G4int columnId, G4float value)
{
  return FillNtupleFColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleDColumn(G4int columnId, G4double value)
{
  return FillNtupleDColumn(fFirstNtupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::FillNtupleIColumn(G4int ntupleId, G4int columnId, 
                                               G4int value)
{
  tools::waxml::ntuple::column<int>* column 
    = GetNtupleIColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::FillNtupleIColumn()",
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
G4bool G4XmlAnalysisManager::FillNtupleFColumn(G4int ntupleId, G4int columnId, 
                                               G4float value)
{
  tools::waxml::ntuple::column<float>* column 
    = GetNtupleFColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::FillNtupleFColumn()",
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
G4bool G4XmlAnalysisManager::FillNtupleDColumn(G4int ntupleId, G4int columnId, 
                                               G4double value)
{
  tools::waxml::ntuple::column<double>* column 
    = GetNtupleDColumn(ntupleId, columnId);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << columnId << " does not exist.";
    G4Exception("G4XmlAnalysisManager::FillNtupleDColumn()",
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
G4bool G4XmlAnalysisManager::AddNtupleRow()
{ 
  return AddNtupleRow(fFirstNtupleId);
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::AddNtupleRow(G4int ntupleId)
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fpVerboseL4->Message("add", "ntuple row", description);
  }  
#endif

  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "AddNtupleRow");

  if ( ! ntupleDescription || ! ntupleDescription->fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4XmlAnalysisManager::AddNtupleRow()",
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
tools::histo::h1d*  G4XmlAnalysisManager::GetH1(G4int id, G4bool warn,
                                                G4bool onlyIfActive) const 
{
  return GetH1InFunction(id, "GetH1", warn, onlyIfActive);
}

//_____________________________________________________________________________
tools::histo::h2d*  G4XmlAnalysisManager::GetH2(G4int id, G4bool warn,
                                                G4bool onlyIfActive) const 
{
  return GetH2InFunction(id, "GetH2", warn, onlyIfActive);
}

//_____________________________________________________________________________
G4int  G4XmlAnalysisManager::GetH1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH1NameIdMap.find(name);
  if ( it ==  fH1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4XmlAnalysisManager::GetH1Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int  G4XmlAnalysisManager::GetH2Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH2NameIdMap.find(name);
  if ( it ==  fH2NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4XmlAnalysisManager::GetH2Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
tools::waxml::ntuple* G4XmlAnalysisManager::GetNtuple() const
{
  return GetNtuple(fFirstNtupleId);
}  

//_____________________________________________________________________________
tools::waxml::ntuple* G4XmlAnalysisManager::GetNtuple(G4int ntupleId) const
{
  G4XmlNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtuple");
    
  return ntupleDescription->fNtuple;
}  

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::GetH1Nbins(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Nbins");
  if ( ! h1d ) return 0;
  
  return h1d->axis().bins();
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH1Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Xmin");
  if ( ! h1d ) return 0;
  
  G4HnInformation* info = GetInformation(kH1, id);
  return info->fXFcn(h1d->axis().lower_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH1Xmax(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Xmax");
  if ( ! h1d ) return 0;
  
  G4HnInformation* info = GetInformation(kH1, id);
  return info->fXFcn(h1d->axis().upper_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH1Width(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1XWidth", true, false);
  if ( ! h1d ) return 0;
  
  G4int nbins = h1d->axis().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH1Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH1, id);
  return ( info->fXFcn(h1d->axis().upper_edge()) 
           - info->fXFcn(h1d->axis().lower_edge()))*info->fXUnit/nbins;
}  

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::GetH2Nxbins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NXbins");
  if ( ! h2d ) return 0;
  
  return h2d->axis_x().bins();
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmin");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fXFcn(h2d->axis_x().lower_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2Xmax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmax");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fXFcn(h2d->axis_x().upper_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2XWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XWidth", true, false);
  if ( ! h2d ) return 0;
  
  G4int nbins = h2d->axis_x().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH2, id);
  return ( info->fXFcn(h2d->axis_x().upper_edge()) 
           - info->fXFcn(h2d->axis_x().lower_edge()))*info->fXUnit/nbins;
}  

//_____________________________________________________________________________
G4int G4XmlAnalysisManager::GetH2Nybins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NYbins");
  if ( ! h2d ) return 0;
  
  return h2d->axis_y().bins();
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymin");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fYFcn(h2d->axis_y().lower_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2Ymax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymax");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fYFcn(h2d->axis_y().upper_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double G4XmlAnalysisManager::GetH2YWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YWidth", true, false);
  if ( ! h2d ) return 0;
  
  G4int nbins = h2d->axis_y().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH2, id);
  return ( info->fYFcn(h2d->axis_y().upper_edge()) 
           - info->fYFcn(h2d->axis_y().lower_edge()))*info->fYUnit/nbins;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH1Title(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1Title");
  if ( ! h1d ) return false;
  
  return h1d->set_title(title);
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH1XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1XAxisTitle");
  if ( ! h1d ) return false;
  
  h1d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH1YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1YAxisTitle");
  if ( ! h1d ) return false;
  
  h1d->add_annotation(tools::histo::key_axis_y_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH2Title(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2Title");
  if ( ! h2d ) return false;
  
  return h2d->set_title(title);
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2XAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2YAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2ZAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_z_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH1Title(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Title");
  if ( ! h1d ) return "";
  
  return h1d->title();
}  

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH1XAxisTitle(G4int id) const 
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1XAxisTitle");
  if ( ! h1d ) return "";
  
  G4String title;
  G4bool result = h1d->annotation(tools::histo::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h1 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH1XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH1YAxisTitle(G4int id) const 
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1YAxisTitle");
  if ( ! h1d ) return "";
  
  G4String title;
  G4bool result = h1d->annotation(tools::histo::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h1 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH1YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH2Title(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Title");
  if ( ! h2d ) return "";
  
  return h2d->title();
}  


//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH2XAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h2 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH2XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
} 

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH2YAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h2 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH2YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String G4XmlAnalysisManager::GetH2ZAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2ZAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_z_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get z_axis title for h2 id = " << id << ").";
    G4Exception("G4XmlAnalysisManager::GetH2ZAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  
