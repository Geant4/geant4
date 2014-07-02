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
#include "G4RootFileManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4RootRNtupleManager.hh"
#include "G4RootRNtupleDescription.hh"
#include "G4AnalysisVerbose.hh"
#include "G4Threading.hh"

#include <tools/rroot/file>
#include <tools/rroot/streamers>
#include <tools/rroot/fac>
#include <tools/rroot/tree>
#include <tools/rroot/ntuple>

#include <iostream>
#include <cstdio>

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
  fNtupleManager = new G4RootRNtupleManager(fState);
  fFileManager = new G4RootFileManager(fState);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
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
G4bool G4RootAnalysisReader::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;
  
  G4bool result = fH1Manager->Reset();
  finalResult = finalResult && result;

  result = fH2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisReader::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;
  G4bool result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "read analysis file", fileName);
#endif

  result = fFileManager->OpenRFile(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read analysis file", fileName, finalResult);
#endif
  
  return finalResult;
}  
  
//_____________________________________________________________________________
G4int G4RootAnalysisReader::GetH1Impl(const G4String& h1Name, 
                                      const G4String& fileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("get", "h1", h1Name);
#endif

  tools::rroot::file* rfile
    = fFileManager->GetRFile(fileName);

  if ( ! rfile ) {
    if ( ! OpenFile(fileName) ) return -1;
    rfile = fFileManager->GetRFile(fileName);
/*
    if ( ! rfile ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Cannot get open file  " << fileName <<"."; 
      G4Exception("G4RootAnalysisReader::GetH1Impl()",
                "Analysis_WR002", JustWarning, description);
      return -1;
    }
*/
  } 
  
  tools::rroot::key* key = rfile->dir().find_key(h1Name);
  if ( ! key ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Key " << h1Name << " for H1 not found in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetH1Impl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(size);
  if ( ! charBuffer ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get data buffer for H1 " << h1Name << " in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetH1Impl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
  }
  
  G4bool verbose = false;
  tools::rroot::buffer buffer(G4cout, rfile->byte_swap(), size, charBuffer, 
                              key->key_length(), verbose);
  tools::histo::h1d* h1 = tools::rroot::TH1D_stream(buffer);
  G4int id = fH1Manager->AddH1(h1Name, h1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("get", "h1", h1Name, id > -1);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::GetH2Impl(const G4String& h2Name, 
                                      const G4String& fileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("get", "h2", h2Name);
#endif

  tools::rroot::file* rfile
    = fFileManager->GetRFile(fileName);

  if ( ! rfile ) {
    if ( ! OpenFile(fileName) ) return -1;
    rfile = fFileManager->GetRFile(fileName);
/*
    if ( ! rfile ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Cannot get open file  " << fileName <<"."; 
      G4Exception("G4RootAnalysisReader::GetH2Impl()",
                "Analysis_WR002", JustWarning, description);
      return -1;
    }
*/
  }
  
  tools::rroot::key* key = rfile->dir().find_key(h2Name);
  if ( ! key ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Key " << h2Name << " for H2 not found in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetH2Impl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(size);
  if ( ! charBuffer ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get data buffer for H2 " << h2Name << " in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetH2Impl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
  }
  
  G4bool verbose = false;
  tools::rroot::buffer buffer(G4cout, rfile->byte_swap(), size, charBuffer, 
                              key->key_length(), verbose);
  tools::histo::h2d* h2 = tools::rroot::TH2D_stream(buffer);
  G4int id = fH2Manager->AddH2(h2Name, h2);  
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("get", "h2", h2Name, id > -1);
#endif
  
  return id;  
  
}  

//_____________________________________________________________________________
G4int G4RootAnalysisReader::GetNtupleImpl(const G4String& ntupleName, 
                                          const G4String& fileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("get", "ntuple", ntupleName);
#endif

  tools::rroot::file* rfile
    = fFileManager->GetRFile(fileName);

  if ( ! rfile ) {
/*
    G4cout << "go to open read file " << fileName << std::endl;
    if ( ! fFileManager->OpenRFile(fileName) ) {
      G4cout << "open read file failed !!" << std::endl;
      // Add exception
      return -1;
    }
*/
    if ( ! OpenFile(fileName) ) return -1;
    rfile = fFileManager->GetRFile(fileName);
  } 
  
  tools::rroot::key* key = rfile->dir().find_key(ntupleName);
  if ( ! key ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Key " << ntupleName << " for Ntuple not found in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetNtupleImpl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(size);
  if ( ! charBuffer ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get data buffer for Ntuple " << ntupleName << " in file " << fileName; 
    G4Exception("G4RootAnalysisReader::GetNtupleImpl()",
                "Analysis_WR002", JustWarning, description);
    return -1;
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
    G4Exception("G4RootAnalysisReader::GetH1Impl()",
                "Analysis_WR002", JustWarning, description);
                
    delete buffer;
    delete tree;    
    return -1;
  }
  
  tools::rroot::ntuple* rntuple 
    = new tools::rroot::ntuple(*tree); //use the flat ntuple API.
  G4RootRNtupleDescription* rntupleDescription
    = new G4RootRNtupleDescription(rntuple, buffer, fac, tree);

  G4int id = fNtupleManager->SetNtuple(rntupleDescription); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("get", "ntuple", ntupleName, id > -1);
#endif

  return id;
}  
