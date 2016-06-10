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
// $Id: G4RootAnalysisReader.cc 74257 2013-10-02 14:24:55Z gcosmo $

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#include "G4RootAnalysisReader.hh"
#include "G4RootRFileManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4RootRNtupleManager.hh"
#include "G4RootRNtupleDescription.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Threading.hh"

#include <tools/rroot/file>
#include <tools/rroot/streamers>
#include <tools/rroot/fac>
#include <tools/rroot/tree>
#include <tools/rroot/ntuple>

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

G4RootAnalysisReader* G4RootAnalysisReader::fgMasterInstance = 0;
G4ThreadLocal G4RootAnalysisReader* G4RootAnalysisReader::fgInstance = 0;

//_____________________________________________________________________________
G4RootAnalysisReader* G4RootAnalysisReader::Instance()
{
  if ( fgInstance == 0 ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4RootAnalysisReader(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4RootAnalysisReader::G4RootAnalysisReader(G4bool isMaster)
 : G4VAnalysisReader("Root", isMaster),
   fH1Manager(0),
   fH2Manager(0),
   fH3Manager(0),
   fP1Manager(0),
   fP2Manager(0),
   fNtupleManager(0),
   fFileManager(0)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4RootAnalysisReader already exists." 
      << "Cannot create another instance.";
    G4Exception("G4RootAnalysisReader::G4RootAnalysisReader()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // Create managers
  fH1Manager = new G4H1ToolsManager(fState);
  fH2Manager = new G4H2ToolsManager(fState);
  fH3Manager = new G4H3ToolsManager(fState);
  fP1Manager = new G4P1ToolsManager(fState);
  fP2Manager = new G4P2ToolsManager(fState);
  fNtupleManager = new G4RootRNtupleManager(fState);
  fFileManager = new G4RootRFileManager(fState);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetH3Manager(fH3Manager);
  SetP1Manager(fP1Manager);
  SetP2Manager(fP2Manager);
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4RootAnalysisReader::~G4RootAnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = 0;
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
tools::rroot::buffer* G4RootAnalysisReader::GetBuffer(
                                               const G4String& fileName,
                                               const G4String& objectName,
                                               const G4String& inFunction)
{
// Get buffer for reading histogram or profile specified by objectNmae
// for a file specified by fileName; 
// open the file if it was not yet open 

  // Histograms and profiles are not saved per thread
  G4bool isPerThread = false;
  
  // Get or open a file
  tools::rroot::file* rfile
    = fFileManager->GetRFile(fileName, isPerThread);
  if ( ! rfile ) {
    if ( ! fFileManager->OpenRFile(fileName, isPerThread) ) return 0;
    rfile = fFileManager->GetRFile(fileName, isPerThread);
  } 
  
  tools::rroot::key* key 
    = ( ! rfile ) ? 0 : rfile->dir().find_key(objectName);
 
  unsigned int size;
  //char* charBuffer 
  //  = ( ! key ) ? 0 : key->get_object_buffer(size);
  char* charBuffer = 0;
  if ( key ) charBuffer = key->get_object_buffer(size);
  
  if ( ! charBuffer ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get " << objectName << " in file " << fileName; 
    G4Exception(inFunction, "Analysis_WR011", JustWarning, description);
    return 0;
  }  

  G4bool verbose = false;
  return new tools::rroot::buffer(G4cout, rfile->byte_swap(), size, charBuffer, 
                                  key->key_length(), verbose);
}

//_____________________________________________________________________________
G4bool G4RootAnalysisReader::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;
  
  G4bool result = fH1Manager->Reset();
  finalResult = finalResult && result;

  result = fH2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fH3Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP1Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
// 
// protected methods
//

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadH1Impl(const G4String& h1Name, 
                                       const G4String& fileName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h1", h1Name);
#endif

  tools::rroot::buffer* buffer 
    = GetBuffer(fileName, h1Name, "ReadH1Impl");
  if ( ! buffer ) return kInvalidId;
  
  tools::histo::h1d* h1 = tools::rroot::TH1D_stream(*buffer);
  delete buffer;
  
  if ( ! h1 ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Streaming " << h1Name << " in file " << fileName << " failed."; 
    G4Exception("G4RootAnalysisReader::ReadH1Impl", 
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }  
  
  G4int id = fH1Manager->AddH1(h1Name, h1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h1", h1Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadH2Impl(const G4String& h2Name, 
                                       const G4String& fileName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h2", h2Name);
#endif

  tools::rroot::buffer* buffer = GetBuffer(fileName, h2Name, "ReadH2Impl");
  if ( ! buffer ) return kInvalidId;
  
  // if h2Name represents H1, then we get !!
  // tools::rroot::buffer::check_byte_count : object of class "TNamed" read too few bytes (603979762 missing).
  // tools::rroot::buffer::check_byte_count : "TNamed" streamer not in sync with data on file, fix streamer.
  // Segmentation fault (core dumped)
  
  tools::histo::h2d* h2 = tools::rroot::TH2D_stream(*buffer);
  delete buffer;
  
  if ( ! h2 ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Streaming " << h2Name << " in file " << fileName << " failed."; 
    G4Exception("G4RootAnalysisReader::ReadH2Impl", 
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }  
  
  G4int id = fH2Manager->AddH2(h2Name, h2);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h2", h2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadH3Impl(const G4String& h3Name, 
                                       const G4String& fileName,
                                       G4bool /*isUserFileName*/)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h3", h3Name);
#endif

  tools::rroot::buffer* buffer = GetBuffer(fileName, h3Name, "ReadH3Impl");
  if ( ! buffer ) return kInvalidId;
  
  tools::histo::h3d* h3 = tools::rroot::TH3D_stream(*buffer);
  delete buffer;
  
  if ( ! h3 ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Streaming " << h3Name << " in file " << fileName << " failed."; 
    G4Exception("G4RootAnalysisReader::ReadH3Impl", 
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }  
  
  G4int id = fH3Manager->AddH3(h3Name, h3);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h3", h3Name, id > kInvalidId);
#endif
  
  return id;  
/* 
  // not yet available
  return kInvalidId;
*/
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadP1Impl(const G4String& p1Name, 
                                       const G4String& fileName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p1", p1Name);
#endif

  tools::rroot::buffer* buffer = GetBuffer(fileName, p1Name, "ReadP1Impl");
  if ( ! buffer ) return kInvalidId;
  
  tools::histo::p1d* p1 = tools::rroot::TProfile_stream(*buffer);
  delete buffer;
  
  if ( ! p1 ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Streaming " << p1Name << " in file " << fileName << " failed."; 
    G4Exception("G4RootAnalysisReader::ReadP1Impl", 
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }  
  
  G4int id = fP1Manager->AddP1(p1Name, p1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p1", p1Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadP2Impl(const G4String& p2Name, 
                                       const G4String& fileName,
                                       G4bool /*isUserFileName*/)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p2", p2Name);
#endif

  tools::rroot::buffer* buffer = GetBuffer(fileName, p2Name, "ReadP2Impl");
  if ( ! buffer ) return kInvalidId;
  
  tools::histo::p2d* p2 = tools::rroot::TProfile2D_stream(*buffer);
  delete buffer;
  
  if ( ! p2 ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Streaming " << p2Name << " in file " << fileName << " failed."; 
    G4Exception("G4RootAnalysisReader::ReadP2Impl", 
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }  
  
  G4int id = fP2Manager->AddP2(p2Name, p2);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p2", p2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::ReadNtupleImpl(const G4String& ntupleName, 
                                           const G4String& fileName,
                                           G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "ntuple", ntupleName);
#endif

  // Ntuples are saved per thread
  // but do not apply the thread suffix if fileName is provided explicitly
  G4bool isPerThread = true;
  if ( isUserFileName ) isPerThread = false; 
  
  // Get or open a file
  tools::rroot::file* rfile
    = fFileManager->GetRFile(fileName, isPerThread);
  if ( ! rfile ) {
    if ( ! fFileManager->OpenRFile(fileName, isPerThread) ) return kInvalidId;
    rfile = fFileManager->GetRFile(fileName, isPerThread);
  } 
  
  tools::rroot::key* key = rfile->dir().find_key(ntupleName);
  if ( ! key ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Key " << ntupleName << " for Ntuple not found in file " << fileName; 
    G4Exception("G4RootAnalysisReader::ReadNtupleImpl()",
                "Analysis_WR011", JustWarning, description);
    return kInvalidId;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(size);
  if ( ! charBuffer ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get data buffer for Ntuple " << ntupleName << " in file " << fileName; 
    G4Exception("G4RootAnalysisReader::ReadNtupleImpl()",
                "Analysis_WR021", JustWarning, description);
    return kInvalidId;
  }
  
  G4bool verbose = false;
  tools::rroot::buffer* buffer
    = new tools::rroot::buffer(G4cout, rfile->byte_swap(), size, charBuffer, 
                               key->key_length(), verbose);
  tools::rroot::fac* fac 
    = new tools::rroot::fac(*rfile);

  tools::rroot::tree* tree 
    = new tools::rroot::tree(*rfile, *fac);
  if ( ! tree->stream(*buffer) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "TTree streaming failed for Ntuple " << ntupleName << " in file " << fileName; 
    G4Exception("G4RootAnalysisReader::ReadNtupleImpl()",
                "Analysis_WR021", JustWarning, description);
                
    delete buffer;
    delete tree;    
    return kInvalidId;
  }
  
  tools::rroot::ntuple* rntuple 
    = new tools::rroot::ntuple(*tree); //use the flat ntuple API.
  G4RootRNtupleDescription* rntupleDescription
    = new G4RootRNtupleDescription(rntuple, buffer, fac, tree);

  G4int id = fNtupleManager->SetNtuple(rntupleDescription); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "ntuple", ntupleName, id > kInvalidId);
#endif

  return id;
}  
