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

using namespace G4Analysis;

namespace {

//_____________________________________________________________________________
void FileManagerWarning(const G4String& fileName,
                        std::string_view className,
                        std::string_view functionName,
                        G4bool hdf5Warn = true)
{
  if ( GetExtension(fileName) == "hdf5" && ( ! hdf5Warn ) ) return;

  Warn("Cannot get file manager for " + fileName,
       className, functionName);
}

}

//_____________________________________________________________________________
G4GenericFileManager::G4GenericFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state)
{}

//
// private methods
//

//_____________________________________________________________________________
void G4GenericFileManager::CreateFileManager(G4AnalysisOutput output)
{
  Message(kVL4, "create", "file manager", GetOutputName(output));

  auto outputId = static_cast<size_t>(output);
  if ( fFileManagers[outputId] ) {
    Warn("The file manager of " + G4Analysis::GetOutputName(output) +
         " type already exists.",
         fkClass, "CreateFileManager");
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
        Warn("Hdf5 type is not available.", fkClass, "CreateFileManager");
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
      Warn(G4Analysis::GetOutputName(output) + " type is not supported.",
        fkClass, "CreateFileManager");
      return;
  }

  // Pass directory names (only if set)
  if ( ! GetHistoDirectoryName().empty() ) {
    fFileManagers[outputId]->SetHistoDirectoryName(GetHistoDirectoryName());
  }
  if ( ! GetNtupleDirectoryName().empty() ) {
    fFileManagers[outputId]->SetNtupleDirectoryName(GetNtupleDirectoryName());
  }

  Message(kVL3, "create", "file manager", GetOutputName(output));
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
  if (extension.size() == 0u) {
    // use the default
    extension = fDefaultFileType;
  }

  auto output = G4Analysis::GetOutput(extension);
  if ( output == G4AnalysisOutput::kNone ) {
    Warn("The file extension " + extension + "is not supported.",
      fkClass, "GetFileManager");
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
    Warn("Default file manager changed "
         "(old: " +fDefaultFileManager->GetFileType() +
         ", new:" + fileManager->GetFileType() + ")",
         fkClass, "OpenFile");
  }
  fDefaultFileManager = fileManager;
  fDefaultFileType = fileManager->GetFileType();

  Message(kVL4, "open", "analysis file", fileName);

  auto result = true;

  // Save the default file name
  // both in the generic file manager and the output specific one
  result &= SetFileName(fileName);
  result &= fDefaultFileManager->SetFileName(fileName);
  result &= fDefaultFileManager->OpenFile(fileName);

  LockDirectoryNames();
  fIsOpenFile = true;

  Message(kVL1, "open", "analysis file", fileName, result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::OpenFiles()
{
// Open all files regeistered with objects

  Message(kVL4, "open", "analysis files");

  auto result = true;

  // process names registered in base file manager
  for ( const auto& fileName : GetFileNames() ) {
    auto fileManager = GetFileManager(fileName);
    if ( ! fileManager ) {
      FileManagerWarning(fileName, fkClass, "OpenFiles", fHdf5Warn);
      continue;
    }

    // filenames for csv need to be updated
    auto newFileName = fileName;
    if (fileManager == fCsvFileManager) {
      newFileName = fileManager->GetHnFileName(fileName, GetCycle());
    }

    result &= fileManager->CreateFile(newFileName);
  }

  Message(kVL3, "open", "analysis files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::WriteFiles()
{
// Finish write for all files registered with objects

  Message(kVL4, "write", "analysis files");

  auto result = true;

  for ( const auto& fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

    Message(kVL4, "write", fileManager->GetFileType(), "files");

    result &= fileManager->WriteFiles();
  }

  Message(kVL3, "write", "analysis files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::CloseFiles()
{
// Close all files regeistered with objects

  Message(kVL4, "close", "analysis files");

  auto result = true;

  for ( const auto& fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

    Message(kVL4, "close", fileManager->GetFileType(), "files");

    result &= fileManager->CloseFiles();
  }

  fIsOpenFile = false;

  Message(kVL3, "close", "analysis files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::DeleteEmptyFiles()
{
// Close all files regeistered with objects

  Message(kVL4, "delete", "empty files");

  auto result = true;

  for ( const auto& fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

    Message(kVL4, "delete", fileManager->GetFileType(), "empty files");

    result &= fileManager->DeleteEmptyFiles();
  }

  Message(kVL3, "delete", "empty files", "", result);

  return result;
}

//_____________________________________________________________________________
void G4GenericFileManager::Clear()
{
// Clear files data

  for ( const auto& fileManager : fFileManagers ) {
    if ( ! fileManager ) continue;

    fileManager->Clear();
  }
  UnlockDirectoryNames();
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::CreateFile(const G4String& fileName)
{
// New prototype, fully implemented in templated base class

  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerWarning(fileName, fkClass, "CreateFile", fHdf5Warn);
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
    FileManagerWarning(fileName, fkClass, "WriteFile", fHdf5Warn);
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
    FileManagerWarning(fileName, fkClass, "CloseFile", fHdf5Warn);
    return false;
  }

  return fileManager->CloseFile(fileName);
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::SetIsEmpty(const G4String& fileName, G4bool isEmpty)
{
  auto fileManager = GetFileManager(fileName);
  if ( ! fileManager ) {
    FileManagerWarning(fileName, fkClass, "SetIsEmpty", fHdf5Warn);
    return false;
  }

  return fileManager->SetIsEmpty(fileName, isEmpty);
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::SetHistoDirectoryName(const G4String& dirName)
{
  auto result = G4VFileManager::SetHistoDirectoryName(dirName);

  for (auto& fileManager : fFileManagers ) {
    if ( fileManager != nullptr ) {
      result &= fileManager->SetHistoDirectoryName(dirName);
    }
  }
  return result;
}

//_____________________________________________________________________________
G4bool G4GenericFileManager::SetNtupleDirectoryName(const G4String& dirName)
{
  auto result = G4VFileManager::SetNtupleDirectoryName(dirName);

  for (auto& fileManager : fFileManagers ) {
    if ( fileManager != nullptr ) {
      result &= fileManager->SetNtupleDirectoryName(dirName);
    }
  }
  return result;
}

//_____________________________________________________________________________
void G4GenericFileManager::SetDefaultFileType(const G4String& value)
{
  // Check if value correspond to a valid file type
  auto output = G4Analysis::GetOutput(value);
  if ( output == G4AnalysisOutput::kNone ) {
    Warn("The file type " + value + "is not supported.\n" +
         "The default type " + fDefaultFileType + " will be used.",
         fkClass, "SetDeafultFileType");
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
      Warn("Failed to create ntuple file manager of " +
           G4Analysis::GetOutputName(output) + " type.\n" + failure,
           fkClass, "CreateNtupleFileManager");
  }

  return vNtupleFileManager;
}
