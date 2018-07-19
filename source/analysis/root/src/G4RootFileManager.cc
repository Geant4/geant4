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
// $Id: G4RootFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4RootFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/zlib"

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootFileManager::G4RootFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fFile(nullptr),
   fHistoDirectory(nullptr),
   fNtupleDirectory(nullptr),
   fNofNtupleFiles(0),
   fNtupleFiles(),
   fMainNtupleDirectories(),
   fBasketSize(0)
{}

//_____________________________________________________________________________
G4RootFileManager::~G4RootFileManager()
{}

//
// private methods
//

//_____________________________________________________________________________
G4bool G4RootFileManager::OpenNtupleFiles()
{
  auto finalResult = true;

  for ( auto i = 0; i < fNofNtupleFiles; i++ ) {
    
    auto name = GetNtupleFileName(i);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "main ntuple file", name);
#endif

    // create new file
    auto rfile = std::make_shared<tools::wroot::file>(G4cout, name);
    rfile->add_ziper('Z', tools::compress_buffer);
    rfile->set_compression(fState.GetCompressionLevel());
    
    if ( ! rfile->is_open() ) {
      G4ExceptionDescription description;
      description << "      " << "Cannot open file " << name;
      G4Exception("G4RootAnalysisManager::OpenFile()",
                  "Analysis_W001", JustWarning, description);
      finalResult = false;
    }

    // Do not create directory if extra ntuple files
    // auto result = CreateNtupleDirectory(rfile);
    // finalResult = finalResult && result;
    tools::wroot::directory* directory = &rfile->dir();
    if ( fNtupleDirectoryName != "" ) {
      directory = rfile->dir().mkdir(fNtupleDirectoryName);
      if ( ! directory ) {
        G4ExceptionDescription description;
        description << "      " 
                    << "cannot create directory " << fNtupleDirectoryName;
        G4Exception("G4RootFileManager::OpenNtupleFiles()",
                    "Analysis_W001", JustWarning, description);
        directory = &fFile->dir();
      }
    }       

    fNtupleFiles.push_back(rfile);
    fMainNtupleDirectories.push_back(directory);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("create", "main ntuple file", name);
#endif

  }

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4RootFileManager::WriteFile(std::shared_ptr<tools::wroot::file> rfile, 
#ifdef G4VERBOSE
                                    const G4String& fileName)
#else
                                    const G4String& /*fileName*/)
#endif
{
  // Do nothing if there is no file
  if ( ! fIsOpenFile ) return true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write", "file", fileName);
#endif

  unsigned int n;
  auto result = rfile->write(n);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("write", "file", fileName, result);
#endif

  return result;
}
  
//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFile(std::shared_ptr<tools::wroot::file> rfile, 
#ifdef G4VERBOSE
                                    const G4String& fileName)
#else
                                    const G4String& /*fileName*/)
#endif
{
  // Do nothing if there is no file
  if ( ! fIsOpenFile ) return true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "file", fileName);
#endif

  rfile->close();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("close", "file", fileName, true);
#endif

  return true;
}
  
//
// public methods
//
//_____________________________________________________________________________
G4bool G4RootFileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;
  auto name = GetFullFileName();
  
  // delete previous file if exists
  //if ( fFile ) delete fFile;

  // create new file
  fFile = std::make_shared<tools::wroot::file>(G4cout, name);
  fFile->add_ziper('Z',tools::compress_buffer);
  fFile->set_compression(fState.GetCompressionLevel());
  
  if ( ! fFile->is_open() ) {
    fFile = nullptr;
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4RootAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }

  // Create directories
  if ( ! CreateHistoDirectory() ) return false;
  if ( ! CreateNtupleDirectory() ) return false;

  // Open ntuple files
  OpenNtupleFiles();

  fLockFileName = true;
  fLockHistoDirectoryName = true;
  fLockNtupleDirectoryName = true;

  fIsOpenFile = true;

  return true;
}  

//_____________________________________________________________________________
G4bool G4RootFileManager:: WriteFile()
{
  auto finalResult = true;

  auto result = WriteFile(fFile, GetFullFileName());
  finalResult = finalResult && result;

  auto counter = 0;
  for ( auto ntupleFile : fNtupleFiles ) {
    result = WriteFile(ntupleFile, GetNtupleFileName(counter++));
    finalResult = finalResult && result;
  }
  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFile()
{
  auto finalResult = true;

  auto result = CloseFile(fFile, GetFullFileName());
  finalResult = finalResult && result;

  auto counter = 0;
  for ( auto ntupleFile : fNtupleFiles ) {
    result = CloseFile(ntupleFile, GetNtupleFileName(counter++));
    finalResult = finalResult && result;
  }

  fLockFileName = false;
  fIsOpenFile = false;

  return finalResult;
} 
   
//_____________________________________________________________________________
G4bool G4RootFileManager::CreateHistoDirectory()
{
  if ( fHistoDirectoryName == "" ) {
    // Do not create a new directory if its name is not set
    fHistoDirectory = &(fFile->dir());
    return true;
  }  
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "directory for histograms", fHistoDirectoryName);
#endif
  
  fHistoDirectory = fFile->dir().mkdir(fHistoDirectoryName);
  if ( ! fHistoDirectory ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot create directory " << fHistoDirectoryName;
    G4Exception("G4RootFileManager::CreateHistoDirectory()",
              "Analysis_W001", JustWarning, description);
    return false;       
  }       
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("create", "directory for histograms", fHistoDirectoryName);
  }    
#endif
  return true;
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CreateNtupleDirectory()
{
  if ( fNtupleDirectoryName == "" ) {
    // Do not create a new directory if its name is not set
    fNtupleDirectory = &(fFile->dir());
    return true;
  }  
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "directory for ntuples", fNtupleDirectoryName);
#endif

  fNtupleDirectory = fFile->dir().mkdir(fNtupleDirectoryName);
  if ( ! fNtupleDirectory ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot create directory " << fNtupleDirectoryName;
    G4Exception("G4RootFileManager::CreateNtupleDirectory()",
                "Analysis_W001", JustWarning, description);
    return false;       
  }       
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("create", "directory for ntuples", fNtupleDirectoryName);
  }    
#endif
  return true;
}

//_____________________________________________________________________________
std::shared_ptr<tools::wroot::file> 
G4RootFileManager::GetNtupleFile(G4int index) const
{ 
  if ( index==0 && ( ! fNtupleFiles.size() ) ) return fFile;

  if ( index < 0 || index >= G4int(fNtupleFiles.size()) ) {
    G4String inFunction = "G4RootFileManager::GetNtupleFile()";
    G4ExceptionDescription description;
    description << "      " << "ntuple file " << index << " does not exist.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return nullptr;         
  }

  return fNtupleFiles[index]; 
}

//_____________________________________________________________________________
tools::wroot::directory*
G4RootFileManager::GetMainNtupleDirectory(G4int index) const
{ 
  if ( index==0 && ( ! fMainNtupleDirectories.size() ) ) return fNtupleDirectory;

  if ( index < 0 || index >= G4int(fMainNtupleDirectories.size()) ) {
    G4String inFunction = "G4RootFileManager::GetMainNtupleDirectory()";
    G4ExceptionDescription description;
    description << "      " << "main ntuple directory " << index << " does not exist.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return nullptr;         
  }

  return fMainNtupleDirectories[index]; 
}

