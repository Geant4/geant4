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
// $Id: G4RootAnalysisManager.cc 103532 2017-04-13 14:00:35Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootAnalysisManager.hh"
#include "G4RootFileManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootPNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <cstdio>

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH2Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH3Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging P1 profiles
  G4Mutex mergeP1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging P2 profiles
  G4Mutex mergeP2Mutex = G4MUTEX_INITIALIZER;
}  

G4RootAnalysisManager* G4RootAnalysisManager::fgMasterInstance = nullptr;
G4ThreadLocal G4RootAnalysisManager* G4RootAnalysisManager::fgInstance = nullptr;

//_____________________________________________________________________________
G4RootAnalysisManager* G4RootAnalysisManager::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4RootAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::IsInstance()
{
  return ( fgInstance != 0 );
}    

//_____________________________________________________________________________
G4RootAnalysisManager::G4RootAnalysisManager(G4bool isMaster)
 : G4ToolsAnalysisManager("Root", isMaster),
   fNofNtupleFiles(0),
   fNtupleMergeMode(G4NtupleMergeMode::kNone),
   fNtupleManager(nullptr),
   fSlaveNtupleManager(nullptr),
   fFileManager(nullptr)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4RootAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4RootAnalysisManager::G4RootAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // File manager
  fFileManager = std::make_shared<G4RootFileManager>(fState);
  SetFileManager(fFileManager);
  fFileManager->SetBasketSize(fgkDefaultBasketSize);

  // Do not merge ntuples by default
  // Merging may require user code migration as analysis manager
  // must be created both on master and workers.
  auto mergeNtuples = false;
  SetNtupleMergingMode(mergeNtuples, fNofNtupleFiles);

  // Create ntuple managers
  CreateNtupleManagers();
}
  
//_____________________________________________________________________________
G4RootAnalysisManager::~G4RootAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootAnalysisManager::SetNtupleMergingMode(G4bool mergeNtuples, 
                                                 G4int nofNtupleFiles)

{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("set", "ntuple merging mode", "");
#endif

  auto canMerge = true;

  // Illegal situations
  if ( mergeNtuples && ( ! G4Threading::IsMultithreadedApplication() ) ) {
    if ( nofNtupleFiles > 0 ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Merging ntuples is not applicable in sequential application." 
        << G4endl 
        << "      " << "Setting was ignored.";
        G4Exception("G4RootAnalysisManager::SetNtupleMergingMode()",
                  "Analysis_W013", JustWarning, description);
    }
    canMerge = false;      
  }

  // Illegal situations
  if ( mergeNtuples && G4Threading::IsMultithreadedApplication() &&
       ( ! fgMasterInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "Merging ntuples requires G4AnalysisManager instance on master." 
      << G4endl 
      << "      " << "Setting was ignored.";
      G4Exception("G4RootAnalysisManager::SetNtupleMergingMode()",
                "Analysis_W013", JustWarning, description);
    canMerge = false;      
  }

  G4String mergingMode;
  if ( ( ! mergeNtuples ) || ( ! canMerge ) ) {
    fNtupleMergeMode = G4NtupleMergeMode::kNone;
    mergingMode = "G4NtupleMergeMode::kNone";      
  }
  else {
    // Set the number of reduced ntuple files
    // G4int nofThreads = G4Threading::GetNumberOfThreads();
    fNofNtupleFiles = nofNtupleFiles;
  
    // Check the number of reduced ntuple files
    // if ( fNofNtupleFiles < 0 || fNofNtupleFiles > nofThreads ) {
    if ( fNofNtupleFiles < 0  ) {
      G4ExceptionDescription description;
      description 
        << "      " << "Number of reduced files must be [0, nofThreads]."
        << G4endl 
        << "      " << "Cannot set  " <<  nofNtupleFiles
        // << " files when nofThreads is " << nofThreads << G4endl   
        << " files" << G4endl   
        << "      " << "Ntuples will be merged in a single file.";
        G4Exception("G4RootAnalysisManager::SetNtupleMergingMode()",
                  "Analysis_W013", JustWarning, description);
      fNofNtupleFiles = 0;
    }
  
    // if ( fNofNtupleFiles == nofThreads ) {
    //   // add warning that no merging will be applied
    //   fNtupleMergeMode = G4NtupleMergeMode::kNone;
    //   fNofNtupleFiles = 0;
    //   mergingMode = "G4NtupleMergeMode::kNone";
    // }
    // else {
    //   G4bool isMaster = ! G4Threading::IsWorkerThread();
    //   if ( isMaster ) {
    //     fNtupleMergeMode = G4NtupleMergeMode::kMain;
    //     mergingMode = "G4NtupleMergeMode::kMain";
    //   } else {
    //     fNtupleMergeMode = G4NtupleMergeMode::kSlave;    
    //     mergingMode = "G4NtupleMergeMode::kSlave";
    //   }
    // }

    // Forced merging mode
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    if ( isMaster ) {
      fNtupleMergeMode = G4NtupleMergeMode::kMain;
      mergingMode = "G4NtupleMergeMode::kMain";
    } else {
      fNtupleMergeMode = G4NtupleMergeMode::kSlave;    
      mergingMode = "G4NtupleMergeMode::kSlave";
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()
      ->Message("set", "ntuple merging mode", mergingMode);
#endif
}

//_____________________________________________________________________________
void G4RootAnalysisManager::ClearNtupleManagers()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("clear", "ntuple managers", "");
#endif

  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave ) {
    // Do not reset master ntuple manager
    delete fNtupleManager;
    fNtupleManager = nullptr;
    // SetNtupleManager(fNtupleManager);
  }

  delete fSlaveNtupleManager;
  fSlaveNtupleManager = nullptr;
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("clear", "ntuple managers", "");
#endif
}

//_____________________________________________________________________________
void G4RootAnalysisManager::CreateNtupleManagers()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple managers", "");
#endif

  switch ( fNtupleMergeMode )
  {
    case G4NtupleMergeMode::kNone:
      fNtupleManager = new G4RootNtupleManager(fState);
      fNtupleManager->SetFileManager(fFileManager);
      SetNtupleManager(fNtupleManager);
      break;

    case G4NtupleMergeMode::kMain: {
      G4int nofMainManagers = fNofNtupleFiles;
      if ( ! nofMainManagers ) nofMainManagers = 1;
             // create one manager if merging required into the histos & profiles files
      fNtupleManager = new G4RootNtupleManager(fState, nofMainManagers);
      fNtupleManager->SetFileManager(fFileManager);
      SetNtupleManager(fNtupleManager);
      break;
    }

    case G4NtupleMergeMode::kSlave:
      fNtupleManager = fgMasterInstance->fNtupleManager;
        // The master class is used only in Get* functions
      auto mainNtupleManager 
        = fNtupleManager->GetMainNtupleManager(GetNtupleFileNumber()); 
      fSlaveNtupleManager = new G4RootPNtupleManager(mainNtupleManager, fState); 
      SetNtupleManager(fSlaveNtupleManager);
      break;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("create", "ntuple managers", "");
#endif
}

//_____________________________________________________________________________
G4int G4RootAnalysisManager::GetNtupleFileNumber()
{
  if ( ! fNofNtupleFiles ) return 0;

  G4int nofMainManagers = fNofNtupleFiles;
  if ( ! nofMainManagers ) nofMainManagers = 1;

  // Debug - check G4Threading::GetNumberOfRunningWorkerThreads()
  G4cout << "In GetNtupleFileNumber: "
         << G4Threading::GetNumberOfRunningWorkerThreads() << G4endl;

  auto fileNumber = G4Threading::G4GetThreadId() % nofMainManagers;
  return fileNumber;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH1()
{
  auto h1Vector = fH1Manager->GetH1Vector();
  auto hnVector = fH1Manager->GetHnVector();

  if ( ! h1Vector.size() ) return true;

  auto result = true;

  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h1Vector, hnVector, directory, "h1");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH1(&mergeH1Mutex);
    fgMasterInstance->fH1Manager->AddH1Vector(h1Vector);
    lH1.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH2()
{
  auto h2Vector = fH2Manager->GetH2Vector();
  auto hnVector = fH2Manager->GetHnVector();

  if ( ! h2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h2Vector, hnVector, directory, "h2");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH2(&mergeH2Mutex);
    fgMasterInstance->fH2Manager->AddH2Vector(h2Vector);
    lH2.unlock();
  }
  
  return result;  
}  

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH3()
{
  auto h3Vector = fH3Manager->GetH3Vector();
  auto hnVector = fH3Manager->GetHnVector();

  if ( ! h3Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h3Vector, hnVector, directory, "h3");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH3(&mergeH3Mutex);
    fgMasterInstance->fH3Manager->AddH3Vector(h3Vector);
    lH3.unlock();
  }
  
  return result;  
}  

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteP1()
{
  auto p1Vector = fP1Manager->GetP1Vector();
  auto hnVector = fP1Manager->GetHnVector();

  if ( ! p1Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(p1Vector, hnVector, directory, "p1");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lP1(&mergeP1Mutex);
    fgMasterInstance->fP1Manager->AddP1Vector(p1Vector);
    lP1.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteP2()
{
  auto p2Vector = fP2Manager->GetP2Vector();
  auto hnVector = fP2Manager->GetHnVector();

  if ( ! p2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(p2Vector, hnVector, directory, "p2");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lP2(&mergeP2Mutex);
    fgMasterInstance->fP2Manager->AddP2Vector(p2Vector);
    lP2.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteNtuple()
{
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone ) return true;
  
  auto finalResult = true;

  G4String ntupleType;
  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) ntupleType = "main ntuples";
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) ntupleType = "slave ntuples";

#ifdef G4VERBOSE 
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("merge", ntupleType, "");
#endif

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Merge();
    finalResult = result && finalResult;
  }  
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) {
    auto result = fSlaveNtupleManager->Merge();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("merge", ntupleType, "");
#endif
  
  return finalResult;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;
  
  auto result = G4ToolsAnalysisManager::Reset();
  finalResult = finalResult && result;
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone || 
       fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    result = fNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  

  finalResult = finalResult && result;
  
  return finalResult;
}

// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  auto finalResult = true;
  auto result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone )  {

#ifdef G4VERBOSE
    G4String name = fFileManager->GetFullFileName();
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "analysis file", name);
#endif

    result = fFileManager->OpenFile(fileName);
    finalResult = finalResult && result;
    
    fNtupleManager->SetNtupleDirectory(fFileManager->GetNtupleDirectory());
    fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "analysis file", name, finalResult);
#endif
  
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {

#ifdef G4VERBOSE
    G4String name = fFileManager->GetFullFileName();
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "main analysis file", name);
#endif

    fFileManager->SetNofNtupleFiles(fNofNtupleFiles);
    result = fFileManager->OpenFile(fileName);
    finalResult = finalResult && result;

    fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "main analysis file", name, finalResult);
#endif  
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave )  {
    // No file is open by Slave manager
    fSlaveNtupleManager->CreateNtuplesFromMain();
  }

  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteImpl() 
{

  auto finalResult = true;

  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) || 
         ( ! fH3Manager->IsEmpty() ) || ( ! fP1Manager->IsEmpty() ) || 
         ( ! fP2Manager->IsEmpty() ) ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "No master G4RootAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram/profile data will not be merged.";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W031", JustWarning, description);
  }
  
  // H1
  auto result = WriteH1();
  finalResult = finalResult && result;

  // H2
  result = WriteH2();
  finalResult = finalResult && result;

  // H3
  result = WriteH3();
  finalResult = finalResult && result;

  // P1
  result = WriteP1();
  finalResult = finalResult && result;

  // P2
  result = WriteP2();
  finalResult = finalResult && result;

  // Ntuples
  result = WriteNtuple();
  finalResult = finalResult && result;

  // File
  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave )  {
    result = fFileManager->WriteFile();
    finalResult = finalResult && result;
  }

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii(fFileManager->GetFileName());
    finalResult = finalResult && result;
  }   

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::CloseFileImpl()
{
  auto finalResult = true;

  G4bool isNtupleManagerEmpty = fNtupleManager->IsEmpty();
    // the ntuple decription vector is cleared on Reset()
    // in kNoMergeAfterOpen ntuple manager mode 

  // reset data
  auto result = Reset();
  if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W021", JustWarning, description);
  } 
  finalResult = finalResult && result;

  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave )  {
    // close file
    fFileManager->CloseFile(); 
  }

  // No files clean-up in sequential mode
  if ( ! G4Threading::IsMultithreadedApplication() )  return finalResult;

  // Delete files if empty in MT mode
  if ( ( fState.GetIsMaster() && 
         fH1Manager->IsEmpty() && fH2Manager->IsEmpty() && fH3Manager->IsEmpty() &&
         fP1Manager->IsEmpty() && fP2Manager->IsEmpty() && isNtupleManagerEmpty ) ||
       ( ( ! fState.GetIsMaster() ) && isNtupleManagerEmpty &&
             fNtupleMergeMode == G4NtupleMergeMode::kNone ) ) {
    result = ! std::remove(fFileManager->GetFullFileName());
    //  std::remove returns 0 when success
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Removing file " 
                  << fFileManager->GetFullFileName() << " failed";
      G4Exception("G4XmlAnalysisManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
    }            
    finalResult = finalResult && result;
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()
        ->Message("delete", "empty file", fFileManager->GetFullFileName());
#endif
  }

  return finalResult;
}

//
// public methods
//

//_____________________________________________________________________________
void G4RootAnalysisManager::SetNtupleMerging(G4bool mergeNtuples, 
                                             G4int nofNtupleFiles,
                                             unsigned int basketSize)

{
  // Keep basketSize in file manager
  fFileManager->SetBasketSize(basketSize);

  // Set ntuple merging mode 
  SetNtupleMergingMode(mergeNtuples, nofNtupleFiles);

  // Clear existing managers
  ClearNtupleManagers();  

  // Re-create managers
  CreateNtupleManagers();
}
