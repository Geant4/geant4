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

#include "tools/wroot/file"
#include "tools/rroot/file"
#include <tools/gzip_buffer>

#include <iostream>
#include <cstdio>

//_____________________________________________________________________________
G4RootFileManager::G4RootFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fFile(0),
   fHistoDirectory(0),
   fProfileDirectory(0),
   fNtupleDirectory(0)
{
}

//_____________________________________________________________________________
G4RootFileManager::~G4RootFileManager()
{  
  delete fFile;
}

//
// public methods
//
//_____________________________________________________________________________
G4bool G4RootFileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;
  G4String name = GetFullFileName();
  
  // delete previous file if exists
  if ( fFile ) delete fFile;

  // create new file
  fFile = new tools::wroot::file(G4cout, name);
  fFile->add_ziper('Z',tools::gzip_buffer);
  fFile->set_compression(9);
  
  if ( ! fFile->is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4RootAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }

  // Create directories
  if ( ! CreateHistoDirectory() ) return false;
  if ( ! CreateProfileDirectory() ) return false;
  if ( ! CreateNtupleDirectory() ) return false;

  fLockFileName = true;
  fLockHistoDirectoryName = true;
  fLockProfileDirectoryName = true;
  fLockNtupleDirectoryName = true;

  return true;
}  
  
//_____________________________________________________________________________
G4bool G4RootFileManager:: WriteFile()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write", "file", GetFullFileName());
#endif

  unsigned int n;
  G4bool result = fFile->write(n);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("write", "file", GetFullFileName(), result);
#endif

  return result;
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFile()
{
  // close file
  fFile->close();  
  fLockFileName = false;

  return true;
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
G4bool G4RootFileManager::CreateProfileDirectory()
{
  if ( fProfileDirectoryName == "" ) {
    // Do not create a new directory if its name is not set
    fProfileDirectory = &(fFile->dir());
    return true;
  }  
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "directory for profiles", fProfileDirectoryName);
#endif
  
  fProfileDirectory = fFile->dir().mkdir(fProfileDirectoryName);
  if ( ! fProfileDirectory ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot create directory " << fProfileDirectoryName;
    G4Exception("G4RootFileManager::CreateProfileDirectory()",
              "Analysis_W001", JustWarning, description);
    return false;       
  }       
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("create", "directory for profiles", fProfileDirectoryName);
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
