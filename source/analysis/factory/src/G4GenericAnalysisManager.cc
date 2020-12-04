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
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4NtupleBookingManager.hh"
#include "G4VNtupleFileManager.hh"

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

// mutex in a file scope

namespace {

//Mutex to lock master manager when merging histograms
G4Mutex mergeHnMutex = G4MUTEX_INITIALIZER;

void WriteHnException(const G4String& hnType, G4int id) 
{
  G4String inFunction = "G4GenericAnalysisManager::Write" + hnType;
  G4ExceptionDescription description;
  description << "Failed to get " << hnType << " id " << id << G4endl;
  G4Exception(inFunction, "Analysis_W022", JustWarning, description);
}

}

G4GenericAnalysisManager* G4GenericAnalysisManager::fgMasterInstance = nullptr;
G4ThreadLocal G4GenericAnalysisManager* G4GenericAnalysisManager::fgInstance = nullptr;

//_____________________________________________________________________________
G4GenericAnalysisManager* G4GenericAnalysisManager::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4GenericAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::IsInstance()
{
  return ( fgInstance != nullptr );
}    

//_____________________________________________________________________________
G4GenericAnalysisManager::G4GenericAnalysisManager(G4bool isMaster)
 : G4ToolsAnalysisManager("", isMaster),
   fFileManager(nullptr),
   fNtupleFileManager(nullptr)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4GenericAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4GenericAnalysisManager::G4GenericAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // File manager
  fFileManager = std::make_shared<G4GenericFileManager>(fState);
  SetFileManager(fFileManager);  
}

//_____________________________________________________________________________
G4GenericAnalysisManager::~G4GenericAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
void G4GenericAnalysisManager::CreateNtupleFileManager(const G4String& fileName)
{
  if ( fNtupleFileManager ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "The ntuple file manager already exists.";
    G4Exception("G4GenericAnalysisManager::CreateNtupleFileManager",
                "Analysis_W002", JustWarning, description);
    return;
  }

  auto fileType = GetExtension(fileName);
  auto output = G4Analysis::GetOutput(fileType);
  if ( output == G4AnalysisOutput::kNone ) {
    G4ExceptionDescription description;
    description 
      << "The file type " << fileType << "is not supported.";
    G4Exception("G4GenericAnalysisManager::CreateNtupleFileManager",
                "Analysis_W051", JustWarning, description);
    return;
  }

  // Set file type to booked ntuples
  fNtupleBookingManager->SetFileType(fileType);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("create", "ntuple file manager", fileType);
  }
#endif

  fNtupleFileManager = fFileManager->CreateNtupleFileManager(output);
  if (fNtupleFileManager) {
    fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
  }

  if ( fNtupleFileManager->IsNtupleMergingSupported() ) {
    // set merginng
    fNtupleFileManager->SetNtupleMerging(fMergeNtuples, fNofNtupleFiles);
    fNtupleFileManager->SetNtupleRowWise(fNtupleRowWise, fNtupleRowMode);
    fNtupleFileManager->SetBasketSize(fBasketSize);
    fNtupleFileManager->SetBasketEntries(fBasketEntries);
  }
  else if ( fIsNtupleMergingSet && fMergeNtuples ) {
    G4ExceptionDescription description;
    description
      << "      " << "Ntuple merging is not available with "
      << fileType << " output." << G4endl
      << "      " << "Setting is ignored.";
    G4Exception("G4GenericAnalysisManager::CreateNtupleFileManager",
      "Analysis_W041", JustWarning, description);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("create", "ntuple file manager", fileType, true);
  }
#endif
}

// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  // Add file name extension, if missing
  auto fullFileName = fileName;
  if ( ! GetExtension(fileName).size() ) {
    // G4cout << "File type is not defined, using default: " << fFileManager->GetDefaultFileType() << G4endl;
    fullFileName = fileName + "." + fFileManager->GetDefaultFileType();
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("open (generic)", "file", fileName);
  }
#endif

  // Create ntuple file manager if there are booked ntuples
  if (! fNtupleFileManager) {
    CreateNtupleFileManager(fullFileName);
  }

  // Create ntuple manager
  // and set it to base class which takes then its ownership
  if (fNtupleFileManager) {
    SetNtupleManager(fNtupleFileManager->CreateNtupleManager());
  }

  auto finalResult = true;

  // Open file for histograms/profiles
  auto result = fFileManager->OpenFile(fullFileName);
  finalResult = finalResult && result;

  // Open ntuple file(s) and create ntuples from bookings
  if (fNtupleFileManager) {
    result = fNtupleFileManager->ActionAtOpenFile(fullFileName);
    finalResult = finalResult && result;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("open (generic)", "file", fileName, finalResult);
  }
#endif
  return finalResult;
}
 
//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::WriteImpl() 
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("write (generic)", "files", "");
  }
#endif

  auto finalResult = true;
  auto result = true;

  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) || 
         ( ! fH3Manager->IsEmpty() ) || ( ! fP1Manager->IsEmpty() ) || 
         ( ! fP2Manager->IsEmpty() ) ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "No master G4GenericAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram/profile data will not be merged.";
      G4Exception("G4GenericAnalysisManager::Write()",
                "Analysis_W031", JustWarning, description);
  }
  
  if ( G4Threading::IsWorkerThread() )  {
    result = Merge();
    finalResult = finalResult && result;
  }
  else {
    // GO ON

    // Open all files registered with objects
    fFileManager->OpenFiles();

    // Write all histograms/profile on master
    //
    result = fFileManager->WriteT(fH1Manager->GetH1Vector(), fH1Manager->GetHnVector());
    finalResult = finalResult && result;

    result = fFileManager->WriteT(fH2Manager->GetH2Vector(), fH2Manager->GetHnVector());
    finalResult = finalResult && result;

    result = fFileManager->WriteT(fH3Manager->GetH3Vector(), fH3Manager->GetHnVector());
    finalResult = finalResult && result;

    result = fFileManager->WriteT(fP1Manager->GetP1Vector(), fP1Manager->GetHnVector());
    finalResult = finalResult && result;

    result = fFileManager->WriteT(fP2Manager->GetP2Vector(), fP2Manager->GetHnVector());
    finalResult = finalResult && result;
  }

  // Ntuples
  if (fNtupleFileManager) {
    result = fNtupleFileManager->ActionAtWrite();
    finalResult = finalResult && result;
  }

  // File
  result = fFileManager->WriteFiles();
  finalResult = finalResult && result;

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii(fFileManager->GetFileName());
    finalResult = finalResult && result;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("write (generic)", "files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::CloseFileImpl(G4bool reset)
{
  auto finalResult = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("close (generic)", "files", "");
  }
#endif

  if (fNtupleFileManager) {
    auto result = fNtupleFileManager->ActionAtCloseFile(reset);
    finalResult = finalResult && result;
  }

  // close file
  auto result = fFileManager->CloseFiles();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Closing files failed";
    G4Exception("G4GenericAnalysisManager::CloseFile()",
              "Analysis_W021", JustWarning, description);
  }
  finalResult = finalResult && result;

  // delete empty files
  result = fFileManager->DeleteEmptyFiles();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Deleting empty files failed";
    G4Exception("G4GenericAnalysisManager::CloseFile()",
              "Analysis_W021", JustWarning, description);
  }
  finalResult = finalResult && result;

  if ( reset ) {
    result = Reset();
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4GenericAnalysisManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
    }
  }
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("close (generic)", "files", "", finalResult);
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;

  auto result = G4ToolsAnalysisManager::Reset();
  finalResult = finalResult && result;

  result = fNtupleFileManager->Reset();
  finalResult = result && finalResult;

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4GenericAnalysisManager::Merge() 
{
  // Nothing to be done on master
  if ( ! G4Threading::IsWorkerThread() ) return false;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("merge (generic) on worker", "histograms", "");
  }
#endif

  // The worker manager just adds its histograms to the master
  fH1Manager->Merge(mergeHnMutex, fgMasterInstance->fH1Manager);
  fH2Manager->Merge(mergeHnMutex, fgMasterInstance->fH2Manager);
  fH3Manager->Merge(mergeHnMutex, fgMasterInstance->fH3Manager);
  fP1Manager->Merge(mergeHnMutex, fgMasterInstance->fP1Manager);
  fP2Manager->Merge(mergeHnMutex, fgMasterInstance->fP2Manager);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    fState.GetVerboseL3()->Message("merge (generic) on worker", "histograms", "", true);
  }
#endif

  return true;
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
    WriteHnException("H1", id);
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
    WriteHnException("H2", id);
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
    WriteHnException("H3", id);
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
    WriteHnException("P1", id);
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
    WriteHnException("P2", id);
    return false;
  } 

  auto p2Name = GetP2Name(id);
  return fFileManager->WriteTExtra<tools::histo::p2d>(fileName, p2d, p2Name);
}
