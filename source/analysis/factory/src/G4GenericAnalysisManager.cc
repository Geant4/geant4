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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4GenericAnalysisManager.hh"
#include "G4GenericFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4NtupleBookingManager.hh"
#include "G4VNtupleFileManager.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Exception.hh"

using namespace G4Analysis;
using std::to_string;

// mutex in a file scope

namespace {

//_____________________________________________________________________________
void WriteHnWarning(const G4String& hnType, G4int id,
                    std::string_view inClass,
                    std::string_view inFunction)
{
  Warn("Failed to get " + hnType + " id " + to_string(id), inClass, inFunction);
}

}

//_____________________________________________________________________________
G4GenericAnalysisManager* G4GenericAnalysisManager::Instance()
{
  static G4ThreadLocalSingleton<G4GenericAnalysisManager> instance;
  fgIsInstance = true;
  return instance.Instance();
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::IsInstance()
{
  return fgIsInstance;
}

//_____________________________________________________________________________
G4GenericAnalysisManager::G4GenericAnalysisManager()
 : G4ToolsAnalysisManager("")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // File manager
  fFileManager = std::make_shared<G4GenericFileManager>(fState);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4GenericAnalysisManager::~G4GenericAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgIsInstance = false;
}

//
// private methods
//

//_____________________________________________________________________________
void G4GenericAnalysisManager::CreateNtupleFileManager(const G4String& fileName)
{
  if ( fNtupleFileManager ) {
    Warn("The ntuple file manager already exists.",
      fkClass, "CreateNtupleFileManager");
    return;
  }

  auto fileType = GetExtension(fileName);
  auto output = G4Analysis::GetOutput(fileType);
  if ( output == G4AnalysisOutput::kNone ) {
    Warn("The file type " + fileType + "is not supported.",
      fkClass, "CreateNtupleFileManager");
    return;
  }

  // Set file type to booked ntuples
  fNtupleBookingManager->SetFileType(fileType);

  Message(kVL4, "create", "ntuple file manager", fileType);

  fNtupleFileManager = fFileManager->CreateNtupleFileManager(output);
  if (fNtupleFileManager) {
    fNtupleFileManager->SetBookingManager(fNtupleBookingManager);

    if ( fNtupleFileManager->IsNtupleMergingSupported() ) {
      // set merginng
      fNtupleFileManager->SetNtupleMerging(fMergeNtuples, fNofNtupleFiles);
      fNtupleFileManager->SetNtupleRowWise(fNtupleRowWise, fNtupleRowMode);
      fNtupleFileManager->SetBasketSize(fBasketSize);
      fNtupleFileManager->SetBasketEntries(fBasketEntries);
    }
    else if ( fIsNtupleMergingSet && fMergeNtuples ) {
      Warn("Ntuple merging is not available with " + fileType + " output.\n" +
           "Setting is ignored.",
           fkClass, "CreateNtupleFileManager");
    }
  }

  Message(kVL3, "create", "ntuple file manager", fileType);
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  Message(kVL4, "open", "file", fileName);

  // Add file name extension, if missing
  auto fullFileName = fileName;
  if ( ! GetExtension(fileName).size() ) {
    auto defaultFileType = fFileManager->GetDefaultFileType();
    // G4cout << "File type is not defined, using default: " << defaultFileType << G4endl;
    if ( ! defaultFileType.size() ) {
      G4Exception("G4GenericAnalysisManager::OpenFileImpl", "Analysis_F001",
        FatalException,
        G4String("Cannot open file \"" + fileName + "\".\n"
          "Please, use a file name with an extension or define the default file type\n"
          "via G4AnalysisManager::SetDefaultFileType()"));
    }

    fullFileName = fileName + "." + fFileManager->GetDefaultFileType();
  }

  // Create ntuple file manager if there are booked ntuples
  if (! fNtupleFileManager) {
    CreateNtupleFileManager(fullFileName);
  }

  // Create ntuple manager
  // and set it to base class which takes then its ownership
  if (fNtupleFileManager) {
    SetNtupleManager(fNtupleFileManager->CreateNtupleManager());
  }

  auto result = true;

  // Open file for histograms/profiles
  result &= fFileManager->OpenFile(fullFileName);

  // Open ntuple file(s) and create ntuples from bookings
  if (fNtupleFileManager) {
    result &= fNtupleFileManager->ActionAtOpenFile(fullFileName);
  }

  Message(kVL3, "open", "file", fileName);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteImpl()
{
  Message(kVL4, "write", "files");

  auto result = true;
  if ( G4Threading::IsWorkerThread() )  {
    result &= G4ToolsAnalysisManager::Merge();
  }
  else {
    // Open all files registered with objects
    fFileManager->OpenFiles();

    // Write all histograms/profile on master
    result &= G4ToolsAnalysisManager::WriteImpl();
  }

  // Ntuples
  if (fNtupleFileManager) {
    result &= fNtupleFileManager->ActionAtWrite();
  }

  // File
  result &= fFileManager->WriteFiles();

  // Write ASCII if activated
  if ( IsAscii() ) {
    result &= WriteAscii(fFileManager->GetFileName());
  }

  Message(kVL3, "write", "files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::CloseFileImpl(G4bool reset)
{
  Message(kVL4, "close", "files");

  auto result = true;
  if (fNtupleFileManager) {
    result &= fNtupleFileManager->ActionAtCloseFile(reset);
  }

  // close file
  if ( ! fFileManager->CloseFiles() ) {
    Warn("Closing files failed", fkClass, "CloseFileImpl");
    result = false;
  }

  // delete empty files
  if ( ! fFileManager->DeleteEmptyFiles() ) {
    Warn("Deleting empty files failed", fkClass, "CloseFileImpl");
    result = false;
  }

  if ( reset ) {
    if ( ! Reset() ) {
      Warn("Resetting data failed", fkClass, "CloseFileImpl");
      result = false;
    }
  }

  Message(kVL3, "close", "files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::ResetImpl()
{
// Reset histograms and ntuple

  Message(kVL4, "reset", "");

  auto result = true;
  result &= G4ToolsAnalysisManager::ResetImpl();
  if ( fNtupleFileManager != nullptr ) {
    result &= fNtupleFileManager->Reset();
  }

  Message(kVL3, "reset", "", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteH1(G4int id, const G4String& fileName)
{
  // Experimental extra write

  // Do not write histo on worker (redundant and fails in hdf5 )
  // If default file is not used, users have to call Merge from their code
  if ( G4Threading::IsWorkerThread() ) return false;

  auto h1d = GetH1(id, false);
  if ( ! h1d ) {
    WriteHnWarning("H1", id, fkClass, "WriteH1");
    return false;
  }

  auto h1Name = GetH1Name(id);
  return fFileManager->WriteTExtra<tools::histo::h1d>(fileName, h1d, h1Name);
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteH2(G4int id, const G4String& fileName)
{
  // Experimental extra write

  // Do not write histo on worker (redundant and fails in hdf5 )
  // If default file is not used, users have to call Merge from their code
  if ( G4Threading::IsWorkerThread() ) return false;

  auto h2d = GetH2(id, false);
  if ( ! h2d ) {
    WriteHnWarning("H2", id, fkClass, "WriteH2");
    return false;
  }

  auto h2Name = GetH2Name(id);
  return fFileManager->WriteTExtra<tools::histo::h2d>(fileName, h2d, h2Name);
}
//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteH3(G4int id, const G4String& fileName)
{
  // Experimental extra write

  // Do not write histo on worker (redundant and fails in hdf5 )
  // If default file is not used, users have to call Merge from their code
  if ( G4Threading::IsWorkerThread() ) return false;

  auto h3d = GetH3(id, false);
  if ( ! h3d ) {
    WriteHnWarning("H3", id, fkClass, "WriteH3");
    return false;
  }

  auto h3Name = GetH3Name(id);
  return fFileManager->WriteTExtra<tools::histo::h3d>(fileName, h3d, h3Name);
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteP1(G4int id, const G4String& fileName)
{
  // Experimental extra write

  // Do not write histo on worker (redundant and fails in hdf5 )
  // If default file is not used, users have to call Merge from their code
  if ( G4Threading::IsWorkerThread() ) return false;

  auto p1d = GetP1(id, false);
  if ( ! p1d ) {
    WriteHnWarning("P1", id, fkClass, "WriteP1");
    return false;
  }

  auto p1Name = GetP1Name(id);
  return fFileManager->WriteTExtra<tools::histo::p1d>(fileName, p1d, p1Name);
}
//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteP2(G4int id, const G4String& fileName)
{
  // Experimental extra write

  // Do not write histo on worker (redundant and fails in hdf5 )
  // If default file is not used, users have to call Merge from their code
  if ( G4Threading::IsWorkerThread() ) return false;

  auto p2d = GetP2(id, false);
  if ( ! p2d ) {
    WriteHnWarning("P2", id, fkClass, "WriteP2");
    return false;
  }

  auto p2Name = GetP2Name(id);
  return fFileManager->WriteTExtra<tools::histo::p2d>(fileName, p2d, p2Name);
}
