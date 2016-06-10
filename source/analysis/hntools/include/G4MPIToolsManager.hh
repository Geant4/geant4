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
// $Id$

// The manager class for MPI applications.

// Author: Ivana Hrivnacova, 25/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4MPIToolsManager_h
#define G4MPIToolsManager_h 1

#include "G4AnalysisManagerState.hh"
#include "G4HnInformation.hh"
#include "G4ios.hh"

#include <tools/impi_world>
#include <tools/histo/hmpi>

#include <vector>

class G4MPIToolsManager
{
  public:
    G4MPIToolsManager(const G4AnalysisManagerState& state,
                      tools::histo::hmpi* hmpi)
    : fState(state), fHmpi(hmpi) {}
    virtual ~G4MPIToolsManager() {}
    
  public:
    // methods
    template <typename T>
    G4bool Merge(const std::vector<T*>& htVector,
                 const std::vector<G4HnInformation*>& hnVector);
  private: 
    // methods
    template <typename T>
    G4bool Send(G4int nofActiveT,
                const std::vector<T*>& htVector,
                const std::vector<G4HnInformation*>& hnVector);

    template <typename T>
    G4bool Receive(G4int nofActiveT,
                const std::vector<T*>& htVector,
                const std::vector<G4HnInformation*>& hnVector);

    // data members
    const G4AnalysisManagerState& fState;
    tools::histo::hmpi*  fHmpi;
};

// inline functions

//_____________________________________________________________________________
template <typename T>
G4bool G4MPIToolsManager::Send(G4int nofActiveT,
                               const std::vector<T*>& htVector,
                               const std::vector<G4HnInformation*>& hnVector)
{
  G4bool finalResult = true;

  // send object to destination rank
  // G4cout << "Begin send for " << nofActiveT << G4endl;
  fHmpi->beg_send(nofActiveT);

  // pack objects
  for ( G4int i=0; i<G4int(htVector.size()); ++i ) {
    // skip sending if activation is enabled and HT is inactivated
    auto info = hnVector[i];
    if ( ( fState.GetIsActivation() && ( ! info->GetActivation() ) ) ) continue; 
    // pack histogram for sending
    // G4cout << "Packed " << i << "th T" << G4endl;
    auto ht = htVector[i];
    auto result = fHmpi->pack(*ht);
    finalResult = result && finalResult;
  }

  //G4cout << "Go to send all " << G4endl;
  if ( ! fHmpi->send(fHmpi->rank()) ) {
    G4ExceptionDescription description;
    description  << "    Rank: " << fHmpi->rank() << " : can't send histos.";
    G4Exception("G4H1ToolsManager::Receieve",
                "Analysis_W031", JustWarning, description);
    return false;
  }

  return finalResult;
}

//_____________________________________________________________________________
template <typename T>
G4bool G4MPIToolsManager::Receive(G4int nofActiveT,
                                  const std::vector<T*>& htVector,
                                  const std::vector<G4HnInformation*>& hnVector)
{
  G4int commSize;
  G4bool result = fHmpi->comm_size(commSize);
  if ( ! result ) {
    G4ExceptionDescription description;
    description 
      << "    Failed to get MPI commander size." << G4endl
      << "    Merging will not be performed.";
    G4Exception("G4H1ToolsManager::Merge",
              "Analysis_W031", JustWarning, description);
    return false;
  }

  // get objects from source ranks
  for (G4int srank = 0; srank < commSize; ++srank) {

    // skip destination rank
    if ( srank == fHmpi->rank() ) continue;

    // get objects from this source rank
    //G4cout << "Go to wait_histos " << rank << G4endl;
    using class_pointer = std::pair<std::string,void*>;
    std::vector<class_pointer> hs;
    if ( ! fHmpi->wait_histos(srank, hs) ) {
      G4ExceptionDescription description;
      description  << "    wait_histos from " << srank << " : failed.";
      G4Exception("G4H1ToolsManager::Receieve",
              "Analysis_W031", JustWarning, description);
      return false;
    }

    // check that we got the right number of objects
    if ( G4int(hs.size()) != nofActiveT ) {
      G4ExceptionDescription description;
      description << "    srank: " << srank << " : got " << hs.size() << " objects, "
                  << "while " << nofActiveT  << " were exepected." << G4endl;
      G4Exception("G4H1ToolsManager::Receieve",
              "Analysis_W031", JustWarning, description);
      return false;
    }

    // merge the objects to destination rank
    G4int counter = 0;
    for ( G4int i=0; i<G4int(htVector.size()); ++i ) {
      // skip sending if activation is enabled and HT is inactivated
      auto info = hnVector[i];
      if ( ( fState.GetIsActivation() && ( ! info->GetActivation() ) ) ) continue; 
      // merge histograms
      auto ht = htVector[i];
      auto newHt = static_cast<T*>(hs[counter++].second);
      ht->add(*newHt);
    }
  }
  return true;
}


//_____________________________________________________________________________
template <typename T>
inline G4bool G4MPIToolsManager::Merge(const std::vector<T*>& htVector,
                                       const std::vector<G4HnInformation*>& hnVector)
{
  if ( ! htVector.size() ) return true;

  // Get number of objects to be sent
  G4int nofActiveT = 0;
  if ( fState.GetIsActivation() ) {
    // only activated histograms will be treated
    for ( G4int i=0; i<G4int(htVector.size()); ++i ) {
      auto activation = hnVector[i]->GetActivation();
      if ( activation ) ++nofActiveT;
    }
  } else {
      nofActiveT = G4int(htVector.size());
  } 

  if ( ! nofActiveT ) return true;

  G4int commRank;
  if ( ! fHmpi->comm_rank(commRank) ) {
    G4ExceptionDescription description;
    description 
      << "    Failed to get MPI commander rank." << G4endl
      << "    Merging will not be performed.";
    G4Exception("G4H1ToolsManager::Merge",
                "Analysis_W031", JustWarning, description);
    return false;
  }

  G4bool finalResult = true;

  if ( commRank != fHmpi->rank() ) {

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) {
      G4ExceptionDescription description;
      description << "on rank "  << commRank
           << " destination rank: " << fHmpi->rank();
      fState.GetVerboseL4()->Message("mpi send", "Hn|Pn", description);
  }  
#endif

    auto result = Send(nofActiveT, htVector, hnVector);

    finalResult = result && finalResult;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) {
      G4ExceptionDescription description;
      description << "on rank " << commRank
         << " destination rank: " << fHmpi->rank();
      fState.GetVerboseL1()->Message("send", "Hn|Pn", description);
    }  
#endif

  } else {

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) {
      G4ExceptionDescription description;
      description << "on rank " << commRank
           << " destination rank: " << fHmpi->rank();
      fState.GetVerboseL4()->Message("mpi wait_histos", "Hn|Pn", description);
    }  
#endif

    auto result = Receive(nofActiveT, htVector, hnVector);

    finalResult = result && finalResult;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) {
      G4ExceptionDescription description;
      description << "on rank "  << commRank
         << " destination rank: " << fHmpi->rank();
      fState.GetVerboseL1()->Message("mpi wait_histos", "Hn|Pn", description);
    }  
#endif
  }
  return finalResult;
}

#endif
