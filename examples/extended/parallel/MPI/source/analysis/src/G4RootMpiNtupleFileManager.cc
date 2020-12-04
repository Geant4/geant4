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

// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#include "G4RootMpiNtupleFileManager.hh"
#include "G4RootMpiNtupleManager.hh"
#include "G4RootMpiPNtupleManager.hh"

#include <tools/impi>

using std::make_shared;

//_____________________________________________________________________________
G4RootMpiNtupleFileManager::G4RootMpiNtupleFileManager(const G4AnalysisManagerState& state)
 : G4RootNtupleFileManager(state),
   fImpi(nullptr),
   fMpiRank(-1),
   fMpiSize(0),
   fMpiSlaveNtupleManager(nullptr),
   fNtupleBooked(false)
{}

//_____________________________________________________________________________
G4RootMpiNtupleFileManager::~G4RootMpiNtupleFileManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootMpiNtupleFileManager::SetMpiNtupleMergingMode(
                               G4int nofNtupleFiles)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    fState.GetVerboseL2()->Message("set", "mpi ntuple merging mode", "");
  }
#endif

  auto canMerge = true;

  // Illegal situations
  if ( fMpiSize < 2 ) {
    G4ExceptionDescription description;
    description 
      << "Merging ntuples is not applicable on a single rank." << G4endl 
      << "Setting was ignored.";
      G4Exception("G4RootMpiNtupleFileManager::SetMpiNtupleMergingMode()",
                  "Analysis_W013", JustWarning, description);
    canMerge = false;      
  }

  G4String mergingMode;
  if ( ! canMerge ) {
    fNtupleMergeMode = G4NtupleMergeMode::kNone;
    mergingMode = "G4NtupleMergeMode::kNone";      
  }
  else {
    // Set the number of reduced ntuple files
    // (multiple output files are not yet supported)
    fNofNtupleFiles = nofNtupleFiles;
  
    // Forced merging mode
    // MPI
    if ( fMpiRank >= fMpiSize ) {
      // the extra worker
      fNtupleMergeMode = G4NtupleMergeMode::kMain;
      mergingMode = "G4NtupleMergeMode::kMain";
    } else {
      // processing worker
      fNtupleMergeMode = G4NtupleMergeMode::kSlave;
      mergingMode = "G4NtupleMergeMode::kSlave";
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) {
    fState.GetVerboseL1()->Message("set", "ntuple merging mode", mergingMode);
  }
#endif
}


// 
// public methods
//

//_____________________________________________________________________________
void G4RootMpiNtupleFileManager::SetMpiNtupleMerging(tools::impi* impi, 
                                             G4int mpiRank, G4int mpiSize,
                                             G4int nofNtupleFiles)
{
  if ( fIsInitialized ) {
      G4ExceptionDescription description;
      description 
        << "Cannot change merging mode." << G4endl
        << "The function must be called before OpenFile().";
      G4Exception("G4RootMpiNtupleFileManager::SetMpiNtupleMerging",
                  "Analysis_W013", JustWarning, description);
      return;
  }

  // Save MPI merging parameters
  fImpi = impi;
  fMpiRank = mpiRank;
  fMpiSize = mpiSize;

  // Set ntuple merging mode 
  SetMpiNtupleMergingMode(nofNtupleFiles);
}

//_____________________________________________________________________________
std::shared_ptr<G4VNtupleManager> G4RootMpiNtupleFileManager::CreateNtupleManager()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("create", "mpi ntuple managers", "");
  }
#endif

  std::shared_ptr<G4VNtupleManager> activeNtupleManager = nullptr;
  switch ( fNtupleMergeMode )
  {
    case G4NtupleMergeMode::kNone:
      fNtupleManager 
        = make_shared<G4RootNtupleManager>(fState, fBookingManager,
                        0, 0, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;

    case G4NtupleMergeMode::kMain: {
      fNtupleManager 
        = make_shared<G4RootMpiNtupleManager>(fState, fBookingManager,
                        fNtupleRowWise, fNtupleRowMode, fImpi, fMpiSize);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;
    }

    case G4NtupleMergeMode::kSlave: {
      auto destinationRank = fMpiSize;
      fMpiSlaveNtupleManager 
        = make_shared<G4RootMpiPNtupleManager>(fState, fImpi, fMpiRank, destinationRank);
      activeNtupleManager = fMpiSlaveNtupleManager;
      break;
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    G4String mergeMode;
    switch ( fNtupleMergeMode ) {
      case G4NtupleMergeMode::kNone:
        mergeMode = "";
        break;
      case G4NtupleMergeMode::kMain:
        mergeMode = "main ";
        break;
      case G4NtupleMergeMode::kSlave:
        mergeMode = "slave ";
        break;
    }
    fState.GetVerboseL3()->Message("create", mergeMode + "mpi ntuple managers", "");
  }
#endif

  fIsInitialized = true;

  return activeNtupleManager;
}

//_____________________________________________________________________________
G4bool G4RootMpiNtupleFileManager::ActionAtOpenFile(const G4String& fileName)
{
  auto finalResult = true;

  // No MPI merging, call base class
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone )  {
    return G4RootNtupleFileManager::ActionAtOpenFile(fileName);
  }

  if ( ! fNtupleBooked ) {

#ifdef G4VERBOSE
    G4String objectType = "analysis file";
    if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) {
      objectType = "main analysis file";
    }
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()->Message("open", objectType, fileName);
    }
#endif

    if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
      // Creating files is triggered from CreateNtuple   
      fNtupleManager->CreateNtuplesFromBooking(
        fBookingManager->GetNtupleBookingVector());
    }

    if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave )  {
      // G4cout << "Slave: Go to create ntuples from booking" << G4endl;
      // No file is open by Slave manager
      fMpiSlaveNtupleManager->CreateNtuplesFromBooking(
        fBookingManager->GetNtupleBookingVector());
    }

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() )  {
      fState.GetVerboseL1()->Message("open", objectType, fileName, finalResult);
    }
#endif

    fNtupleBooked = true;
  }

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4RootMpiNtupleFileManager::ActionAtWrite()
{
  // No MPI merging, call base class
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone  ) {
    return G4RootNtupleFileManager::ActionAtWrite();
  }
  
  auto finalResult = true;

  G4String ntupleType;
  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) ntupleType = "main ntuples";
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) ntupleType = "slave ntuples";

#ifdef G4VERBOSE 
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("merge", ntupleType, "");
  }
#endif

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Merge();
    finalResult = result && finalResult;
  }  
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) {
    auto result = fMpiSlaveNtupleManager->Merge();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) {
    fState.GetVerboseL1()->Message("merge", ntupleType, "");
  }
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootMpiNtupleFileManager::ActionAtCloseFile(G4bool reset)
{
  // No MPI merging, call base class
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone )  {
    return G4RootNtupleFileManager::ActionAtCloseFile(reset);
  }

  auto finalResult = true;

  // reset data
  if ( reset ) {
    auto result = Reset();
    if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "Resetting data failed";
        G4Exception("G4RootNtupleFileManager::Write()",
                  "Analysis_W021", JustWarning, description);
    } 
    finalResult = finalResult && result;
  }

  // close files
  if ( fNtupleMergeMode != G4NtupleMergeMode::kSlave )  {
    auto result = CloseNtupleFiles();
    finalResult = finalResult && result;
  }

  // MT not yet supported - no files clean-up
  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootMpiNtupleFileManager::Reset()
{
// Reset ntuples

  // No MPI merging, call base class
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone  ) {
    return G4RootNtupleFileManager::Reset();
  }
  
  auto finalResult = true;
  
  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  

  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave )  {
    auto result = fMpiSlaveNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  

  return finalResult;
}
