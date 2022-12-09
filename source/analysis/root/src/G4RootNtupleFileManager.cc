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

// Author: Ivana Hrivnacova, 15/09/2020 (ivana@ipno.in2p3.fr)

#include "G4RootNtupleFileManager.hh"
#include "G4RootFileManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootPNtupleManager.hh"
#include "G4VAnalysisManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"

using namespace G4Analysis;
using std::make_shared;
using std::to_string;



G4RootNtupleFileManager* G4RootNtupleFileManager::fgMasterInstance = nullptr;

//_____________________________________________________________________________
G4RootNtupleFileManager::G4RootNtupleFileManager(const G4AnalysisManagerState& state)
 : G4VNtupleFileManager(state, "root")
{
  if ( G4Threading::IsMasterThread() ) fgMasterInstance = this;

  // Do not merge ntuples by default
  // Merging may require user code migration as analysis manager
  // must be created both on master and workers.
  auto mergeNtuples = false;
  SetNtupleMergingMode(mergeNtuples, fNofNtupleFiles);
}

//_____________________________________________________________________________
G4RootNtupleFileManager::~G4RootNtupleFileManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
}

//
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleMergingMode(G4bool mergeNtuples,
                                                 G4int nofNtupleFiles)

{
  Message(kVL4, "set", "ntuple merging mode");

  auto canMerge = true;

  // Illegal situations
  if ( mergeNtuples && ( ! G4Threading::IsMultithreadedApplication() ) ) {
    Warn("Merging ntuples is not applicable in sequential application.\n"
         "Setting was ignored.",
         fkClass, "SetNtupleMergingMode");
    canMerge = false;
  }

  // Illegal situations
  if (mergeNtuples && G4Threading::IsMultithreadedApplication() && (fgMasterInstance == nullptr)) {
    Warn("Merging ntuples requires G4AnalysisManager instance on master.\n"
         "Setting was ignored.",
         fkClass, "SetNtupleMergingMode");
    canMerge = false;
  }

  G4String mergingMode;
  if ( ( ! mergeNtuples ) || ( ! canMerge ) ) {
    fNtupleMergeMode = G4NtupleMergeMode::kNone;
    mergingMode = "G4NtupleMergeMode::kNone";
  }
  else {
    // Set the number of reduced ntuple files
    fNofNtupleFiles = nofNtupleFiles;

    // Check the number of reduced ntuple files
    if ( fNofNtupleFiles < 0  ) {
      Warn("Number of reduced files must be [0, nofThreads].\n"
           "Cannot set  " + to_string(nofNtupleFiles) +" files.\n" +
           "Setting was ignored.",
           fkClass, "SetNtupleMergingMode");
      fNofNtupleFiles = 0;
    }

    // Forced merging mode
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    // G4bool isMaster = fState.GetIsMaster();
    if ( isMaster ) {
      fNtupleMergeMode = G4NtupleMergeMode::kMain;
      mergingMode = "G4NtupleMergeMode::kMain";
    } else {
      fNtupleMergeMode = G4NtupleMergeMode::kSlave;
      mergingMode = "G4NtupleMergeMode::kSlave";
    }
  }

  Message(kVL2, "set", "ntuple merging mode", mergingMode);
}

//_____________________________________________________________________________
G4int G4RootNtupleFileManager::GetNtupleFileNumber()
{
  if (fNofNtupleFiles == 0) return 0;

  G4int nofMainManagers = fNofNtupleFiles;
  if (nofMainManagers == 0) nofMainManagers = 1;

  auto fileNumber = G4Threading::G4GetThreadId() % nofMainManagers;
  return fileNumber;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::CloseNtupleFiles()
{
  // Take into account main ntuple files if present
  auto mainNumber = ( fNofNtupleFiles > 0 ) ? 0 : -1;

  auto result = true;
  auto ntupleVector = fNtupleManager->GetNtupleDescriptionVector();
  for ( auto ntupleDescription : ntupleVector) {
    for (G4int i = mainNumber; i < fNofNtupleFiles; ++i ) {
      result &= fFileManager->CloseNtupleFile(ntupleDescription, i);
    }
  }
  return result;
}

//
// public methods
//

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleMerging(G4bool mergeNtuples,
                                               G4int  nofNtupleFiles)

{
  if ( fIsInitialized ) {
    Warn("Cannot change merging mode.\n"
         "The function must be called before OpenFile().",
         fkClass, "SetNtupleMerging");
    return;
  }

  // Set ntuple merging mode
  SetNtupleMergingMode(mergeNtuples, nofNtupleFiles);
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetNtupleRowWise(G4bool rowWise, G4bool rowMode)
{

  // Print info even when setting makes no effect
  // (as we do not get printed the default setting in the output)
  G4String rowWiseMode;
  if ( rowWise ) {
    rowWiseMode = "row-wise with extra branch";
  }
  else if ( rowMode ) {
    rowWiseMode = "row-wise";
  }
  else {
    rowWiseMode = "column-wise";
  }

  Message(kVL1, "set", "ntuple merging row mode", rowWiseMode);

  // Do nothing if the mode is not changed
  if ( fNtupleRowWise == rowWise && fNtupleRowMode == rowMode ) return;

  fNtupleRowWise = rowWise;
  fNtupleRowMode = rowMode;

  if ( fNtupleManager ) {
    fNtupleManager->SetNtupleRowWise(rowWise, rowMode);
  }

  if ( fSlaveNtupleManager ) {
    fSlaveNtupleManager->SetNtupleRowWise(rowWise, rowMode);
  }
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetBasketSize(unsigned int basketSize)
{
  fFileManager->SetBasketSize(basketSize);
}

//_____________________________________________________________________________
void G4RootNtupleFileManager::SetBasketEntries(unsigned int basketEntries)
{
  fFileManager->SetBasketEntries(basketEntries);
}

//_____________________________________________________________________________
std::shared_ptr<G4VNtupleManager> G4RootNtupleFileManager::CreateNtupleManager()
{
// Create and return ntuple manager.
// If ntuple merging is activated, managers of different types are created
// on master/worker.

  Message(kVL4, "create", "ntuple manager");

  // Check that file manager and anaysis manager are set !

  std::shared_ptr<G4VNtupleManager> activeNtupleManager = nullptr;
  switch ( fNtupleMergeMode )
  {
    case G4NtupleMergeMode::kNone:
      fNtupleManager
        = make_shared<G4RootNtupleManager>(
            fState, fBookingManager, 0, 0, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;

    case G4NtupleMergeMode::kMain: {
      G4int nofMainManagers = fNofNtupleFiles;
      if (nofMainManagers == 0) {
        // create one manager if merging required into the histos & profiles files
        nofMainManagers = 1;
      }
      fNtupleManager
        = make_shared<G4RootNtupleManager>(
            fState, fBookingManager, nofMainManagers, fNofNtupleFiles, fNtupleRowWise, fNtupleRowMode);
      fNtupleManager->SetFileManager(fFileManager);
      activeNtupleManager = fNtupleManager;
      break;
    }

    case G4NtupleMergeMode::kSlave:
      fNtupleManager = fgMasterInstance->fNtupleManager;
        // The master class is used only in Get* functions
      auto mainNtupleManager
        = fNtupleManager->GetMainNtupleManager(GetNtupleFileNumber());
      fSlaveNtupleManager
        = make_shared<G4RootPNtupleManager>(
            fState, fBookingManager, mainNtupleManager, fNtupleRowWise, fNtupleRowMode);
      activeNtupleManager = fSlaveNtupleManager;
      break;
  }

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
  Message(kVL3, "create", mergeMode + "ntuple manager");

  fIsInitialized = true;

  return activeNtupleManager;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtOpenFile([[maybe_unused]] const G4String& fileName)
{
  // Check if fNtupleBookingManager is set

  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone ||
       fNtupleMergeMode == G4NtupleMergeMode::kMain )  {

    G4String objectType = "analysis file";
    if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) {
      objectType = "main analysis file";
    }
    Message(kVL4, "open", objectType, fileName);

    // Creating files is triggered from CreateNtuple
    fNtupleManager->CreateNtuplesFromBooking(
      fBookingManager->GetNtupleBookingVector());

    Message(kVL1, "open", objectType, fileName);
  }

  // Creating ntuples from main is triggered from the first Fill call
  // if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave )  {
  //   // G4cout << "Going to create slave ntuples from main" << G4endl;
  //   fSlaveNtupleManager->CreateNtuplesFromMain();
  // }

  return true;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtWrite()
{
  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone ) {
    return true;
  }

  auto result = true;

  G4String ntupleType;
  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain ) ntupleType = "main ntuples";
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) ntupleType = "slave ntuples";

  Message(kVL4, "merge", ntupleType);

  if ( fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    result &= fNtupleManager->Merge();
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) {
    result &= fSlaveNtupleManager->Merge();
  }

  Message(kVL1, "merge", ntupleType, "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::ActionAtCloseFile()
{
  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave)  {
    fSlaveNtupleManager->SetNewCycle(false);
    return true;
  }

  return CloseNtupleFiles();
}

//_____________________________________________________________________________
G4bool G4RootNtupleFileManager::Reset()
{
// Reset ntuples

  auto result = true;

  if ( fNtupleMergeMode == G4NtupleMergeMode::kNone ||
       fNtupleMergeMode == G4NtupleMergeMode::kMain )  {
    result &= fNtupleManager->Reset();
  }

  if ( fNtupleMergeMode == G4NtupleMergeMode::kSlave ) {
    fSlaveNtupleManager->Reset();
  }

  return result;
}
