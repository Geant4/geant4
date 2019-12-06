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

#include "G4RootMpiAnalysisManager.hh"
#include "G4RootMpiNtupleManager.hh"
#include "G4RootMpiPNtupleManager.hh"

#include <tools/impi>

//_____________________________________________________________________________
G4RootMpiAnalysisManager::G4RootMpiAnalysisManager(G4bool isMaster)
 : G4RootAnalysisManager(isMaster),
   fMpiNtupleMergeMode(G4MpiNtupleMergeMode::kNone),
   fMpiSlaveNtupleManager(nullptr)
{}

//_____________________________________________________________________________
G4RootMpiAnalysisManager::~G4RootMpiAnalysisManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootMpiAnalysisManager::CreateMpiNtupleManagers(
                              tools::impi* impi, G4int mpiRank, G4int mpiSize)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "mpi ntuple managers", "");
#endif

  switch ( fMpiNtupleMergeMode )
  {
    case G4MpiNtupleMergeMode::kNone:
      fNtupleManager 
        = new G4RootNtupleManager(fState, 0, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      SetNtupleManager(fNtupleManager);
      break;

    case G4MpiNtupleMergeMode::kMain: {
      fNtupleManager 
        = new G4RootMpiNtupleManager(fState, fNtupleRowWise, fNtupleRowMode,
                                     impi, mpiSize);
      fNtupleManager->SetFileManager(fFileManager);
      SetNtupleManager(fNtupleManager);
      break;
    }

    case G4MpiNtupleMergeMode::kSlave: {
      auto destinationRank = mpiSize;
      fMpiSlaveNtupleManager 
        = new G4RootMpiPNtupleManager(fState, impi, mpiRank, destinationRank);
      SetNtupleManager(fMpiSlaveNtupleManager);
      break;
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) 
    fState.GetVerboseL3()->Message("create", "mpi ntuple managers", "");
#endif
}

//_____________________________________________________________________________
void G4RootMpiAnalysisManager::SetMpiNtupleMergingMode(
                               G4int mpiRank, G4int mpiSize,
                               G4int nofNtupleFiles)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("set", "mpi ntuple merging mode", "");
#endif

  auto canMerge = true;

  // Illegal situations
  if ( mpiSize < 2 ) {
    G4ExceptionDescription description;
    description 
      << "      " << "Merging ntuples is not applicable on a single rank." 
      << G4endl 
      << "      " << "Setting was ignored.";
      G4Exception("G4RootMpiAnalysisManager::SetMpiNtupleMergingMode()",
                "Analysis_W013", JustWarning, description);
    canMerge = false;      
  }

  G4String mergingMode;
  if ( ! canMerge ) {
    fMpiNtupleMergeMode = G4MpiNtupleMergeMode::kNone;
    mergingMode = "G4MpiNtupleMergeMode::kNone";      
  }
  else {
    // Set the number of reduced ntuple files
    // (multiple output files are not yet supported)
    fNofNtupleFiles = nofNtupleFiles;
  
    // Forced merging mode
    // MPI
    if ( mpiRank >= mpiSize ) {
      // the extra worker
      fMpiNtupleMergeMode = G4MpiNtupleMergeMode::kMain;
      mergingMode = "G4MpiNtupleMergeMode::kMain";
    } else {
      // processing worker
      fMpiNtupleMergeMode = G4MpiNtupleMergeMode::kSlave;
      mergingMode = "G4MpiNtupleMergeMode::kSlave";
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()
      ->Message("set", "ntuple merging mode", mergingMode);
#endif
}


//_____________________________________________________________________________
void G4RootMpiAnalysisManager::SetMpiNtupleMerging(tools::impi* impi, 
                                             G4int mpiRank, G4int mpiSize,
                                             G4int nofNtupleFiles)
{
  // G4cout << "SetMpiNtupleMerging: "
  //        << impi << ", "
  //        << mpiRank << ","
  //        << mpiSize << ","
  //        << nofNtupleFiles << G4endl;

  // fImpi = impi;

  // Set ntuple merging mode 
  SetMpiNtupleMergingMode(mpiRank, mpiSize, nofNtupleFiles);

  // Clear existing managers
  ClearNtupleManagers();  

  // Re-create managers
  CreateMpiNtupleManagers(impi, mpiRank, mpiSize);
}

// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  // No MPI merging, call base class
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kNone )  {
    return G4RootAnalysisManager::OpenFileImpl(fileName);
  }

  auto finalResult = true;
  auto result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kMain )  {

#ifdef G4VERBOSE
    G4String name = fFileManager->GetFullFileName();
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "main ntuple file", name);
#endif

    fFileManager->SetNofNtupleFiles(fNofNtupleFiles);
    result = fFileManager->OpenFile(fileName);
    finalResult = finalResult && result;

    fNtupleManager->SetNtupleDirectory(fFileManager->GetNtupleDirectory());

    G4cout << "Main: Go to create ntuples from booking " << G4endl;
    fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "main ntuple file", name, finalResult);
#endif  
  }

  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kSlave )  {

#ifdef G4VERBOSE
    G4String name = fFileManager->GetFullFileName();
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("open", "file", name);
#endif
    result = fFileManager->OpenFile(fileName);
    finalResult = finalResult && result;

    // fNtupleManager->SetNtupleDirectory(fFileManager->GetNtupleDirectory());

    G4cout << "Slave: Go to create ntuples from booking" << G4endl;
    // No file is open by Slave manager
    fMpiSlaveNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("open", "file", name, finalResult);
#endif  
  }

  return finalResult;
}  


//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::CloseFileImpl(G4bool reset)
{
  // No MPI merging, call base class
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kNone )  {
    return G4RootAnalysisManager::CloseFileImpl(reset);
  }

  auto finalResult = true;

  // reset data
  if ( reset ) {
    auto result = Reset();
    if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "Resetting data failed";
        G4Exception("G4RootAnalysisManager::Write()",
                  "Analysis_W021", JustWarning, description);
    } 
    finalResult = finalResult && result;
  }

  // close file
  fFileManager->CloseFile(); 

  // MT not yet supported - no files clean-up
  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::WriteNtuple()
{
  // No MPI merging, call base class
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kNone  ) {
    return G4RootAnalysisManager::WriteNtuple();
  }
  
  auto finalResult = true;

  G4String ntupleType;
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kMain ) ntupleType = "main ntuples";
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kSlave ) ntupleType = "slave ntuples";

#ifdef G4VERBOSE 
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("merge", ntupleType, "");
#endif

  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kMain )  {
    auto result = fNtupleManager->Merge();
    finalResult = result && finalResult;
  }  
  
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kSlave ) {
    auto result = fMpiSlaveNtupleManager->Merge();
    finalResult = result && finalResult;
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("merge", ntupleType, "");
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::Reset()
{
// Reset histograms and ntuple

  // No MPI merging, call base class
  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kNone  ) {
    return G4RootAnalysisManager::Reset();
  }
  
  auto finalResult = true;
  
  auto result = G4ToolsAnalysisManager::Reset();
  finalResult = finalResult && result;

  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kMain )  {
    result = fNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  

  if ( fMpiNtupleMergeMode == G4MpiNtupleMergeMode::kSlave )  {
    result = fMpiSlaveNtupleManager->Reset(false);
    finalResult = result && finalResult;
  }  
  
  return finalResult;
}
