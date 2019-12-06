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

#include "G4RootMpiNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/mpi_ntuple_row_wise"
#include "tools/wroot/mpi_ntuple_column_wise"

using namespace G4Analysis;

const int kTAG_NTUPLE = 1004;   // This constant is defined in G4MPImanager
                                // (should be passed from the application)

//_____________________________________________________________________________
G4RootMpiNtupleManager::G4RootMpiNtupleManager(
                          const G4AnalysisManagerState& state, 
                          G4bool rowWise, G4bool rowMode,
                          tools::impi* impi, G4int mpiSize)
 : G4RootNtupleManager(state, 0, rowWise, rowMode),
   fImpi(impi),
   fSlaveRanks(),
   fMainRank(0)
{
  for ( G4int rank = 0; rank < mpiSize; rank++ ) {
    fSlaveRanks.push_back(rank);
  }
  fMainRank = mpiSize;
}

//_____________________________________________________________________________
G4RootMpiNtupleManager::~G4RootMpiNtupleManager()
{}

// 
// private methods
//


//_____________________________________________________________________________
G4bool G4RootMpiNtupleManager::Send(G4int id, tools::wroot::ntuple* ntuple)
{
// Pack and send the main ntuple data to the slave ranks

  // G4cout << "Going to send main ntuple data " << G4endl;

  // Get basket sizes    
  std::vector<tools::wroot::branch*> mainBranches;
  ntuple->get_branches(mainBranches);
  std::vector<tools::uint32> basketSizes;
  tools_vforcit(tools::wroot::branch*, mainBranches, it) {
    basketSizes.push_back((*it)->basket_size());
  }

  auto ntupleFile = fFileManager->GetNtupleFile(id);
  tools::uint32 basketSize = fFileManager->GetBasketSize();
  unsigned int basketEntries = fFileManager->GetBasketEntries();

  for ( auto slaveRank : fSlaveRanks ) {
    G4cout << "Going to send main ntuple data to slave rank " << slaveRank << G4endl;

    fImpi->pack_reset();
    if ( ! fImpi->pack(id)) {
      G4cerr << "pack(id) failed." << G4endl; 
      return false;
    }
  
    if ( ! fImpi->bpack(fRowWise) ) {
      G4cerr << "bpack(fRowWise) failed."<< G4endl;
      return false;
    }

    if ( ! fImpi->bpack(ntupleFile->byte_swap())) {
      G4cerr << "bpack(byte_swap) failed." << G4endl; 
      return false;
    }
  
    if ( ! fImpi->pack(ntupleFile->compression()) ) {
      G4cerr << "pack(compression) failed." << G4endl;
      return false; 
    }
  
    if ( ! fImpi->pack(ntupleFile->dir().seek_directory())) {
      // Should be used fNtupleDirectory ??
      G4cerr << "pack(seek) failed." << G4endl;
      return false;
    }

    if ( fRowWise ) {
      if ( ! fImpi->pack(basketSize) ) {
        G4cerr << "pack(basketSize) failed." << G4endl;
        return false;
      }
    } else {
      if ( ! fImpi->bpack(fRowMode) ) {
        G4cerr << "bpack(fRowMode) failed." << G4endl;
        return false;
      }      
      if ( ! fImpi->vpack(basketSizes) ) {
        G4cerr << "vpack(basketSizes) failed." << G4endl;
        return false;
      }
      if ( ! fImpi->pack(basketEntries) ) {
        G4cerr << "pack(basketEntries) failed." << G4endl;
        return false;
      }      
    }

    if ( ! fImpi->send_buffer(slaveRank, kTAG_NTUPLE )) {
      G4cerr << "send_buffer() failed." << G4endl;
      return false;
    }
    fImpi->pack_reset();

    G4cout << "Sent ntuple description to slave on rank " << slaveRank << G4endl;
  }
 
  // G4cout << "Done: send main ntuple data " << G4endl;
  return true;
}

//_____________________________________________________________________________
G4bool G4RootMpiNtupleManager::InitializeRanks()
{
  // G4cout << "G4RootMpiNtupleManager::InitializeRanks" << G4endl;

  auto finalResult = true;

  auto counter = 0;
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {

    // Do not create ntuple if it is inactivated 
    if ( fState.GetIsActivation() && ( ! ntupleDescription->fActivation ) ) continue;

    auto result = Send(counter++, ntupleDescription->fNtuple);
    finalResult = finalResult && result;
  }

  return finalResult;
} 

//_____________________________________________________________________________
G4bool G4RootMpiNtupleManager::WaitBuffer()
{
// Receive the pntuple data from the slave ranks
// For the time being only one ntuple 

  // G4cout << "G4RootMpiNtupleManager::WaitBuffer" << G4endl;

  unsigned long numberOfEndFill = 0;
  
  G4bool verbose = ( fState.GetVerboseL2() );

  while ( true ) { 
    fImpi->pack_reset();

    // loop until receiving end_fill from all ranks
    // G4cout << "G4RootMpiNtupleManager::WaitBuffer entering loop" << G4endl;
    int probe_src;
    if ( ! fImpi->wait_buffer(fMainRank, kTAG_NTUPLE, probe_src, verbose)) {
      G4cerr << "!!! wait_buffer() failed." << std::endl;
      return EXIT_FAILURE;
    }
  
    tools::uint32 protocol;
    if ( ! fImpi->unpack(protocol)) {
      G4cerr << "unpack(protocol) failed."<< G4endl;
      return false;
    }
  
    if ( protocol == tools::wroot::mpi_protocol_basket() ) {

      // G4cout << "G4RootMpiNtupleManager::WaitBuffer got protocol_basket" << G4endl;
  
      // get ntuple Id
      tools::uint32 ntupleId;
      if ( ! fImpi->unpack(ntupleId) ) {
        G4cerr << "unpack(ntuple_id) failed."<< std::endl;
        return false;
      }

      if ( ntupleId >= fNtupleVector.size() ) {
        std::cerr << "!!! unknown ntupleId " << ntupleId << std::endl;
        return false;
      }

      // Main ntuple
      auto mainNtuple = fNtupleVector[ntupleId];

      // add basket to main ntuple
      if ( ! mainNtuple->mpi_add_basket(*fImpi)) {
          std::cerr << "mainNtuple->mpi_add_basket() failed." << std::endl;    
          return EXIT_FAILURE;
      }
    } 
    else if ( protocol == tools::wroot::mpi_protocol_baskets() ) {
      //column_wise and row_mode only.

      // G4cout << "G4RootMpiNtupleManager::WaitBuffer got protocol_baskets" << G4endl;
  
      // get ntuple Id
      tools::uint32 ntupleId;
      if ( ! fImpi->unpack(ntupleId) ) {
        G4cerr << "unpack(ntuple_id) failed."<< std::endl;
        return false;
      }

      if ( ntupleId >= fNtupleVector.size() ) {
        std::cerr << "!!! unknown ntupleId " << ntupleId << std::endl;
        return false;
      }

      // Main ntuple
      auto mainNtuple = fNtupleVector[ntupleId];

      // add basket to main ntuple
      if ( ! mainNtuple->mpi_add_baskets(*fImpi)) {
          std::cerr << "mainNtuple->mpi_add_baskets() failed." << std::endl;    
          return EXIT_FAILURE;
      }
    } 
    else if ( protocol==tools::wroot::mpi_protocol_end_fill() ) {
  
      // G4cout << "G4RootMpiNtupleManager::WaitBuffer got protocol_end_fill" << G4endl;
      
      // get ntuple Id
      tools::uint32 ntupleId;
      if ( ! fImpi->unpack(ntupleId) ) {
        G4cerr << "unpack(ntuple_id) failed."<< std::endl;
        return false;
      }

      if ( ntupleId >= fNtupleVector.size() ) {
        std::cerr << "!!! unknown ntupleId " << ntupleId << std::endl;
        return false;
      }
  
      // Main ntuple
      auto mainNtuple = fNtupleVector[ntupleId];

      // end_fill in main ntuple
      if ( ! mainNtuple->mpi_end_fill(*fImpi) ) {
        G4cerr << "main_ntuple->mpi_end_fill() failed." << std::endl;    
        return false;
      }

      numberOfEndFill++;

      if ( numberOfEndFill == fSlaveRanks.size() ) break;
    } 
    else {
      G4cerr << "unknown protocol " << protocol << G4endl;
      return false;
    }
  }

  return true;
}


// 
// public methods
//

//_____________________________________________________________________________
void G4RootMpiNtupleManager::CreateNtuplesFromBooking()
{
  // Base class actions
  // G4cout << "Going to call CreateNtuplesFromBooking from base class" << G4endl;
  G4TNtupleManager<tools::wroot::ntuple>::CreateNtuplesFromBooking();

  // Initialize ranks
  if ( ! InitializeRanks() ) {
    G4cerr << "InitializeRanks failed." << G4endl;
  }

  // Go to wait buffer mode
  if ( ! WaitBuffer() ) {
    G4cerr << "WaitBuffer failed." << G4endl;
  }
}

//_____________________________________________________________________________
G4bool G4RootMpiNtupleManager::Merge()
{
  G4cout << "G4RootMpiNtupleManager::Merge()" << G4endl;

  auto finalResult = true;

  for ( auto ntupleDescription : fNtupleDescriptionVector ) {

    // Do not create ntuple if it is inactivated 
    if ( fState.GetIsActivation() && ( ! ntupleDescription->fActivation ) ) continue;

    // G4cout << "Go to call merge_number_of_entries" << G4endl;
    ntupleDescription->fNtuple->merge_number_of_entries();
  }

  return finalResult;
}
