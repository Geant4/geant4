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
#include "G4RootMpiNtupleFileManager.hh"
#include "G4RootMpiNtupleManager.hh"
#include "G4RootMpiPNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include <tools/impi>

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootMpiAnalysisManager* G4RootMpiAnalysisManager::Instance()
{
  return fgInstance;
}

//_____________________________________________________________________________
G4RootMpiAnalysisManager::G4RootMpiAnalysisManager(G4bool /*isMaster*/)
 : G4RootAnalysisManager()
{
  fgInstance = this;

  // Reset the ntuple file manager
  fNtupleFileManager.reset();

  // Ntuple file manager
  fNtupleFileManager = std::make_shared<G4RootMpiNtupleFileManager>(fState);
  SetNtupleFileManager(fNtupleFileManager);
  fNtupleFileManager->SetFileManager(fFileManager);
  fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
}

//_____________________________________________________________________________
G4RootMpiAnalysisManager::~G4RootMpiAnalysisManager()
{
  fgInstance = 0;
}

//
// public methods
//

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

  std::static_pointer_cast<G4RootMpiNtupleFileManager>(fNtupleFileManager)
    ->SetMpiNtupleMerging(impi, mpiRank, mpiSize, nofNtupleFiles);
}

// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  if ( fNtupleFileManager->GetMergeMode() == G4NtupleMergeMode::kNone )  {
    return G4RootAnalysisManager::OpenFileImpl(fileName);
  }

  // Create ntuple manager(s)
  // and set it to base class which takes then their ownership
  if ( ! fVNtupleManager ) {
    SetNtupleManager(fNtupleFileManager->CreateNtupleManager());
  }

  auto result = true;

  // Open file
  // In difference from base class a file is open also on slave ranks
  result &= fFileManager->OpenFile(fileName);

  // Open ntuple file(s) and create ntuples from bookings
  result &= fNtupleFileManager->ActionAtOpenFile(fFileManager->GetFullFileName());

  return result;
}  

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::WriteImpl() 
{

  auto result = true;

  // Call base class method
  result &= G4ToolsAnalysisManager::WriteImpl();

  // Write file also on Slave
  // (skipped in base class)
  if ( fNtupleFileManager->GetMergeMode() == G4NtupleMergeMode::kSlave )  {
    // write all open files
    result &= fFileManager->WriteFiles();
  }

  Message(kVL2, "write", "slave files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4RootMpiAnalysisManager::CloseFileImpl(G4bool reset)
{
  auto result = true;

  // Call base class method
  result &= G4RootAnalysisManager::CloseFileImpl(reset);

  // Close file also on Slave
  // (skipped in base class)
  if ( fNtupleFileManager->GetMergeMode() == G4NtupleMergeMode::kSlave )  {
    // close all open files
    result &= fFileManager->CloseFiles();
  }

  Message(kVL2, "close", "slave files", "", result);

  return result;

}
