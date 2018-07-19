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

#include "G4Hdf5FileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/hdf5/h2file"

using namespace G4Analysis;

//_____________________________________________________________________________
const G4String G4Hdf5FileManager::fgkDefaultDirectoryName = "default";

//_____________________________________________________________________________
G4Hdf5FileManager::G4Hdf5FileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fFile(kInvalidId),
   fHistoDirectory(kInvalidId),
   fNtupleDirectory(kInvalidId),
   fBasketSize(0)          // TO DO: check default value !! (433 in test)
{}

//_____________________________________________________________________________
G4Hdf5FileManager::~G4Hdf5FileManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::CreateDirectory(const G4String& directoryType, 
  const G4String& directoryName, hid_t& directory)
{
// Method for both histograms and ntuples directories.
// For histos: CreateDirectory("histograms", fHistoDirectoryName, fHistoDirectory)
// For ntples: CreateDirectory("ntuples", fNtupleDirectoryName, fNtupleDirectory)

  auto newDirectoryName = directoryName;

  if ( newDirectoryName == "" ) {
    // if ( fDefaultDirectory > 0  ) {
    //   // Return the default directory if the name is not set and the default directory
    //   // already exists
    //   directory = fDefaultDirectory;
    //   return true;
    // } else {
      // Create the default directory if the name is not set and the default directory
      // does not yet exist
      newDirectoryName = fgkDefaultDirectoryName;
      newDirectoryName += "_";
      newDirectoryName += directoryType;
  }
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4String message = "directory for ";
    message += directoryType;
    fState.GetVerboseL4()->Message("create", message, newDirectoryName);
  }
#endif
  
  directory = tools_H5Gcreate(fFile, newDirectoryName, 0);
       // 0 seems to be an optional parameter. The web doc does not say what should
       // be the default value but 0 is what is found in examples, and in the code, if we pass 0, clearly some
       // default value is taken.

  if ( directory < 0 ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot create directory " << directoryName;
    G4Exception("G4Hdf5FileManager::CreateDirectory()",
                "Analysis_W001", JustWarning, description);
    return false;       
  }  
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) {
    G4String message = "directory for ";
    message += directoryType;
      fState.GetVerboseL2()->Message("create", message, newDirectoryName);
    }
  }    
#endif
  return true;
}

//_____________________________________________________________________________
#ifdef G4VERBOSE
G4bool G4Hdf5FileManager::WriteDirectory(const G4String& directoryType,
   const G4String& directoryName, hid_t& directory)
#else
G4bool G4Hdf5FileManager::WriteDirectory(const G4String& /*directoryType*/,
   const G4String& directoryName, hid_t& directory)
#endif  
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4String message = "directory for ";
    message += directoryType;
    fState.GetVerboseL4()
      ->Message("write", message, directoryName);
  }
#endif
  
  auto result 
    = tools::hdf5::write_atb(directory, "type", "directory");

  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " 
                << "cannot write directory " << directoryName;
    G4Exception("G4Hdf5FileManager::WriteDirectory()",
                "Analysis_W001", JustWarning, description);
    // close 
    // ::H5Gclose(histos);
    // ::H5Fclose(file);

    return false;
  }  
#ifdef G4VERBOSE
  else {
    if ( fState.GetVerboseL2() ) {
    G4String message = "directory for ";
    message += directoryType;
      fState.GetVerboseL2()->Message("write", message, fHistoDirectoryName);
    }
  }    
#endif

  return result;
}
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;
  auto name = GetFullFileName();
  
  // delete previous file if exists
  //if ( fFile ) delete fFile;

  // create new file
  fFile = ::H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  
  if ( fFile < 0 ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4Hdf5AnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }

  // Create directories
  if ( ! CreateDirectory("histograms", fHistoDirectoryName, fHistoDirectory) ) return false;
  if ( ! CreateDirectory("ntuples", fNtupleDirectoryName, fNtupleDirectory) ) return false;

  // // Open ntuple files
  // OpenNtupleFiles();

  // Write directories
  if ( ! WriteHistoDirectory() ) return false;
  if ( ! WriteNtupleDirectory() ) return false;

  fLockFileName = true;
  fLockHistoDirectoryName = true;
  fLockNtupleDirectoryName = true;

  fIsOpenFile = true;

  return true;
}
 
//_____________________________________________________________________________
G4bool G4Hdf5FileManager::WriteFile() 
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::CloseFile()
{
  // Do nothing if there is no file
  if ( fFile < 0 ) return true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "file", GetFullFileName());
#endif

  if ( fHistoDirectory >= 0 ) {
    ::H5Gclose(fHistoDirectory);
  }
  if ( fNtupleDirectory >= 0 ) {
    ::H5Gclose(fNtupleDirectory);
  }
  ::H5Fclose(fFile);

  fLockFileName = false;
  fIsOpenFile = false;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("close", "file", GetFullFileName(), true);
#endif

  // CHECK
  // the result not returned from hdf5 

  return true;
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::WriteHistoDirectory()
{
  return WriteDirectory("histograms", fHistoDirectoryName, fHistoDirectory);
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::WriteNtupleDirectory()
{
  return WriteDirectory("ntuples", fNtupleDirectoryName, fNtupleDirectory);
}

//_____________________________________________________________________________
void G4Hdf5FileManager::CloseAfterHnWrite()
{
  if ( fHistoDirectory >= 0 ) {
    ::H5Gclose(fHistoDirectory);
  }
  ::H5Fclose(fFile);
}



