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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4GenericFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4CsvFileManager.hh"
#include "G4CsvNtupleFileManager.hh"
#ifdef TOOLS_USE_HDF5
#include "G4Hdf5FileManager.hh"
#include "G4Hdf5NtupleFileManager.hh"
#endif
#include "G4RootFileManager.hh"
#include "G4RootNtupleFileManager.hh"
#include "G4XmlFileManager.hh"
#include "G4XmlNtupleFileManager.hh"

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

namespace {

void FileManagerException(const G4String& fileName, const G4String& functionName,
                          const G4String& exceptionClassification, G4bool hdf5Warn = true)
{
  if ( GetExtension(fileName) == "hdf5" && ( ! hdf5Warn ) ) return;

  G4String where = "G4GenericFileManager::" + functionName;
  G4String what = "Analysis_" + exceptionClassification;
  G4ExceptionDescription description;
  description << "Cannot get file manager for " << fileName;
  G4Exception(where, what, JustWarning, description);
}

}

// static data
const G4String G4GenericFileManager::fgkDefaultFileType = "root";

//_____________________________________________________________________________
G4GenericFileManager::G4GenericFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fDefaultFileType(fgkDefaultFileType),
   fDefaultFileManager(nullptr),
   fFileManagers
     { nullptr, // Csv
       nullptr, // Hdf5
       nullptr, // Generic
       nullptr  // Xml
     },
   fCsvFileManager(nullptr),
#ifdef TOOLS_USE_HDF5
   fHdf5FileManager(nullptr),
#endif
   fRootFileManager(nullptr),
   fXmlFileManager(nullptr),
   fHdf5Warn(true)
{}

//_____________________________________________________________________________
G4GenericFileManager::~G4GenericFileManager()
{}

//
// private methods
//

//_____________________________________________________________________________
void G4GenericFileManager::CreateFileManager(G4AnalysisOutput output)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("create", "file manager", GetOutputName(output));
  }
#endif

  auto outputId = static_cast<size_t>(output);  
  if ( fFileManagers[outputId] ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "The file manager of " << G4Analysis::GetOutputName(output) << " type already exists.";
    G4Exception("G4GenericFileManager::CreateFileManager",
                "Analysis_W002", JustWarning, description);
    return;
  }

  // Create the manager
  switch ( output ) {
    case G4AnalysisOutput::kCsv:
      fCsvFileManager = std::make_shared<G4CsvFileManager>(fState);
      fFileManagers[outputId] = fCsvFileManager;
      break;
    case G4AnalysisOutput::kHdf5: 
#ifdef TOOLS_USE_HDF5
      fHdf5FileManager = std::make_shared<G4Hdf5FileManager>(fState);
      fFileManagers[outputId] = fHdf5FileManager;
#else
      if ( fHdf5Warn) {
        G4ExceptionDescription description;
        description << "Hdf5 type is not available.";
        G4Exception("G4GenericFileManager::CreateFileManager",
                    "Analysis_W051", JustWarning, description);
        fHdf5Warn = false;
      }
#endif
      break;
    case G4AnalysisOutput::kRoot:
      fRootFileManager = std::make_shared<G4RootFileManager>(fState);
      fFileManagers[outputId] = fRootFileManager;
      break;
    case G4AnalysisOutput::kXml:
      fXmlFileManager = std::make_shared<G4XmlFileManager>(fState);
      fFileManagers[outputId] = fXmlFileManager ;
      break;
    case G4AnalysisOutput::kNone:
      G4ExceptionDescription description;
      description 
        << G4Analysis::GetOutputName(output) << " type is not supported.";
      G4Exception("G4GenericFileManager::CreateFileManager",
                  "Analysis_W051", JustWarning, description);
      break;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("create", "file manager", GetOutputName(output));
  }
#endif
}

//_____________________________________________________________________________
std::shared_ptr<G4VFileManager> 
G4GenericFileManager::GetFileManager(G4AnalysisOutput output) const
{
  return fFileManagers[static_cast<size_t>(output)];
}

//_____________________________________________________________________________
std::shared_ptr<G4VFileManager> 
G4GenericFileManager::GetFileManager(const G4String& fileName)
{
  // Get file extension
  G4String extension = GetExtension(fileName);
  if ( ! extension.size() ) {
    // use the default
    extension = fDefaultFileType;
  }

  auto output = G4Analysis::GetOutput(extension);
  if ( output == G4AnalysisOutput::kNone ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "The file extension " << extension << "is not supported.";
    G4Exception("G4GenericFileManager::GetFileManager",
                "Analysis_W051", JustWarning, description);
    return nullptr;
  }

  std::shared_ptr<G4VFileManager> fileManager = GetFileManager(output);
  if ( ! GetFileManager(output) ) {
    CreateFileManager(output);
    fileManager = GetFileManager(output);
  }

  return GetFileManager(output);
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4GenericFileManager::OpenFile(const G4String& fileName)
{
  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) return false;

  if ( fDefaultFileManager && (fDefaultFileManager != fileManager) ) {
    // Print warning if default output changed
    // (maybe be not needed?)
    G4ExceptionDescription description;
    description
      << "Default file manager changed (old: " 
      << fDefaultFileManager->GetFileType()
      << ", new:" << fileManager->GetFileType() << ")";
    G4Exception("G4GenericFileManager::OpenFile",
                "Analysis_W001", JustWarning, description);
  } 
  fDefaultFileManager = fileManager;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("open", "analysis file", fileName);
  }
#endif

  auto finalResult = true;
  auto result = true;

  // Save the default file name 
  // both in the generic file manager and the output specific one
  result = SetFileName(fileName);
  finalResult = finalResult && result;
  result = fDefaultFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

  result = fDefaultFileManager->OpenFile(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) {
    fState.GetVerboseL1()->Message("open", "analysis file", fileName, finalResult);
  }
#endif

  fLockDirectoryNames = true;
  fIsOpenFile = true;

  return finalResult;  
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::OpenFiles()
{
// Open all files regeistered with objects

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("open", "analysis files", "");
  }
#endif

  auto finalResult = true;
  auto result = true;

  // process names registered in base file manager
  for ( auto fileName : GetFileNames() ) {
    auto fileManager = GetFileManager(fileName);
    if ( ! fileManager ) {
      FileManagerException(fileName, "OpenFiles", "W001", fHdf5Warn);
      continue;
    }

    result = fileManager->CreateFile(fileName);
    finalResult = result && finalResult;    
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("open", "analysis files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::WriteFiles()
{
// Finish write for all files regeistered with objects

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("write", "files", "");
  }
#endif

  auto finalResult = true;
  auto result = true;

  for ( auto fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()->Message("write", fileManager->GetFileType(), "files");
    }
#endif

    result = fileManager->WriteFiles();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("write", "files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::CloseFiles()
{
// Close all files regeistered with objects

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("close", "files", "");
  }
#endif

  auto finalResult = true;
  auto result = true;

  for ( auto fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()->Message("close", fileManager->GetFileType(), "files");
    }
#endif

    result = fileManager->CloseFiles();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("close", "files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::DeleteEmptyFiles()
{
// Close all files regeistered with objects

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("delete", "empty files", "");
  }
#endif

  auto finalResult = true;
  auto result = true;

  for ( auto fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()->Message("delete", fileManager->GetFileType(), "files");
    }
#endif

    result = fileManager->DeleteEmptyFiles();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("delete", "empty files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::CreateFile(const G4String& fileName)
{
// New prototype, fully implemented in templated base class

  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerException(fileName, "CreateFile", "W001", fHdf5Warn);
    return false;
  }

  return fileManager->CreateFile(fileName);
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::WriteFile(const G4String& fileName)
{
// New prototype, fully implemented in templated base class

  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerException(fileName, "WriteFile", "W021", fHdf5Warn);
    return false;
  }

  return fileManager->WriteFile(fileName);
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::CloseFile(const G4String& fileName)
{
// New prototype, fully implemented in templated base class

  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerException(fileName, "CloseFile",  "W021", fHdf5Warn);
    return false;
  }

  return fileManager->CloseFile(fileName);
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::SetIsEmpty(const G4String& fileName, G4bool isEmpty)
{
  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerException(fileName, "SetIsEmpty", "W021", fHdf5Warn);
    return false;
  }

  return fileManager->SetIsEmpty(fileName, isEmpty);
}

//_____________________________________________________________________________
void G4GenericFileManager::SetDefaultFileType(const G4String& value)
{
  // Check if value correspond to a valid file type
  auto output = G4Analysis::GetOutput(value);
  if ( output == G4AnalysisOutput::kNone ) {
    G4ExceptionDescription description;
    description
      << "The file type " << value << "is not supported." << G4endl
      << "The default type " << fDefaultFileType << " will be used.";
    G4Exception("G4GenericFileManager::SetDeafultFileType",
                "Analysis_W051", JustWarning, description);
    return;
  }

  fDefaultFileType = value;
}

//_____________________________________________________________________________
std::shared_ptr<G4VNtupleFileManager> 
G4GenericFileManager::CreateNtupleFileManager(G4AnalysisOutput output)
{
  if ( ! GetFileManager(output) ) {
    CreateFileManager(output);
  }

  std::shared_ptr<G4VNtupleFileManager> vNtupleFileManager = nullptr;
  G4String failure;

  switch ( output ) {
    case G4AnalysisOutput::kCsv: {
      auto ntupleFileManager = std::make_shared<G4CsvNtupleFileManager>(fState);
      ntupleFileManager->SetFileManager(fCsvFileManager);
      vNtupleFileManager = ntupleFileManager;
      break;
    } 
    case G4AnalysisOutput::kHdf5: {
#ifdef TOOLS_USE_HDF5
      auto ntupleFileManager = std::make_shared<G4Hdf5NtupleFileManager>(fState);
      ntupleFileManager->SetFileManager(fHdf5FileManager);
      vNtupleFileManager = ntupleFileManager;
#else
      failure = " Hdf5 is not available";
#endif
      break;
    }
    case G4AnalysisOutput::kRoot: {
      auto ntupleFileManager = std::make_shared<G4RootNtupleFileManager>(fState);
      ntupleFileManager->SetFileManager(fRootFileManager);
      vNtupleFileManager = ntupleFileManager;
      break;
    }
    case G4AnalysisOutput::kXml: {
      auto ntupleFileManager = std::make_shared<G4XmlNtupleFileManager>(fState);
      ntupleFileManager->SetFileManager(fXmlFileManager);
      vNtupleFileManager = ntupleFileManager;
      break;
    }
    case G4AnalysisOutput::kNone:
      break;
  }

  if ( ! vNtupleFileManager ) {
      G4ExceptionDescription description;
      description 
        << "      " 
        << "Failed to create ntuple file manager of " << G4Analysis::GetOutputName(output) << " type."
        << failure;
      G4Exception("G4GenericFileManager::CreateNtupleFileManager",
                  "Analysis_W002", JustWarning, description);
  }

  return vNtupleFileManager;
}
