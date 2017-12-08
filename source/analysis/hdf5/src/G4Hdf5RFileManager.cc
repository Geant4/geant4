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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4Hdf5RFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/hdf5/h2file"
#include <tools/zlib>

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

//_____________________________________________________________________________
const G4String G4Hdf5RFileManager::fgkDefaultDirectoryName = "default";

//_____________________________________________________________________________
G4Hdf5RFileManager::G4Hdf5RFileManager(const G4AnalysisManagerState& state)
 : G4BaseFileManager(state),
   fRFiles()
{
}

//_____________________________________________________________________________
G4Hdf5RFileManager::~G4Hdf5RFileManager()
{  
}

// 
// private methods
//

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::OpenRFile(const G4String& fileName,
                                     G4bool isPerThread)
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "read analysis file", name);
#endif

  // create new file
  hid_t newFile = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( newFile < 0 ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << name;
    G4Exception("G4Hdf5RFileManager::OpenFile()",
                "Analysis_WR001", JustWarning, description);
    return false;
  }

  // newFile->add_unziper('Z',tools::decompress_buffer);
  
  // add file in a map
  std::map<G4String, hid_t>::iterator it
    = fRFiles.find(name);
  if ( it != fRFiles.end() ) { 
    it->second = newFile;
  }
  else {
    fRFiles[name] = newFile;
  }    

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read analysis file", name);
#endif

  return true;
}  
  
//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::OpenDirectory(hid_t file, const G4String& directoryName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("open", "read directory", directoryName);
  }
#endif
  
  auto directory = tools_H5Gopen(file, directoryName);
  if ( directory < 0 ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot open directory " << directoryName;
    G4Exception("G4Hdf5RFileManager::OpenDirectory()",
                "Analysis_W001", JustWarning, description);
    return kInvalidId;       
  }  
  else {
#ifdef G4VERBOSE
    if ( fState.GetVerboseL2() ) {
      fState.GetVerboseL2()->Message("open", "read directory", directoryName);
    }
#endif
    return directory;
  }    
}

//_____________________________________________________________________________
hid_t  G4Hdf5RFileManager::GetRDirectory(const G4String& directoryType,
                                        const G4String& fileName, 
                                        const G4String& dirName,
                                        G4bool isPerThread)
{
  // Get or open a file
  auto rfile = GetRFile(fileName, isPerThread);
  if ( rfile < 0 ) {
    // Try to open it if not found in the map 
    if ( ! OpenRFile(fileName, isPerThread) ) return kInvalidId;
    rfile = GetRFile(fileName, isPerThread);
  }

  // Use default directory name if not specified
  auto newDirName = dirName;
  if ( newDirName == "" ) {
    // if ( fDefaultDirectory > 0  ) {
    //   // Return the default directory if the name is not set and the default directory
    //   // already exists
    //   directory = fDefaultDirectory;
    //   return true;
    // } else {
      // Create the default directory if the name is not set and the default directory
      // does not yet exist
      newDirName = fgkDefaultDirectoryName;
      newDirName += "_";
      newDirName += directoryType;
  }

  // Open directory
  return OpenDirectory(rfile, newDirName);
}


//
// public methods
//

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::GetRFile(const G4String& fileName,
                                   G4bool isPerThread) const
{ 
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  std::map<G4String, hid_t>::const_iterator it
    = fRFiles.find(name);
  if  ( it != fRFiles.end() )
    return it->second;
  else {
    return kInvalidId;
  }     
}

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::GetHistoRDirectory(const G4String& fileName, const G4String& dirName,
                                             G4bool isPerThread)
{
  return GetRDirectory("histograms", fileName, dirName, isPerThread);
}

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::GetNtupleRDirectory(const G4String& fileName, const G4String& dirName,
                                              G4bool isPerThread)
{
  return GetRDirectory("ntuples", fileName, dirName, isPerThread);
}


