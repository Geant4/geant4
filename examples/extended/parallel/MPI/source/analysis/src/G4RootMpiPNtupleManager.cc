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

#include "G4RootMpiPNtupleManager.hh"
#include "G4RootFileManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/mpi_ntuple_row_wise"
#include "tools/wroot/mpi_ntuple_column_wise"

const int kTAG_NTUPLE = 1004;   // This constant is defined in G4MPImanager
                                // (should be passed from the application)

//_____________________________________________________________________________
G4RootMpiPNtupleManager::G4RootMpiPNtupleManager(
                           const G4AnalysisManagerState& state, 
                           tools::impi* impi, G4int mpiRank, G4int destinationRank)
 : G4BaseNtupleManager(state),
   fNtupleVector(),
   fImpi(impi),
   fMpiRank(mpiRank),
   fDestinationRank(destinationRank)
{
}

//_____________________________________________________________________________
G4RootMpiPNtupleManager::~G4RootMpiPNtupleManager()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    delete ntupleDescription;
  }   
}

//
// private functions
//

//_____________________________________________________________________________
G4RootMpiPNtupleDescription*  
G4RootMpiPNtupleManager::GetNtupleDescriptionInFunction(
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootMpiPNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return nullptr;         
  }
  
  return fNtupleDescriptionVector[index];
}

//_____________________________________________________________________________
tools::wroot::base_pntuple*  G4RootMpiPNtupleManager::GetNtupleInFunction(
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto ntupleDescription = GetNtupleDescriptionInFunction(id, functionName);
  if ( ! ntupleDescription ) return nullptr;

  if ( ! ntupleDescription->fBasePNtuple ) {
    if ( warn ) {
      G4String inFunction = "G4RootMpiPNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      "  << "ntupleId " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return nullptr;
  }  
  return ntupleDescription->fBasePNtuple;
}

//
// protected functions
//

//_____________________________________________________________________________
void G4RootMpiPNtupleManager::CreateNtuple(G4RootMpiPNtupleDescription* ntupleDescription)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "mpi pntuple", ntupleDescription->fNtupleBooking.name());
#endif

  // Wait for the ntuple data from main

  G4bool verbose = ( fState.GetVerboseL2() );

  G4cout << "Go to wait_buffer from " << fDestinationRank << G4endl;
  fImpi->pack_reset();
  int probe_src;
  if ( ! fImpi->wait_buffer(fMpiRank, fDestinationRank, kTAG_NTUPLE, probe_src, verbose)) {
    G4cerr << "G4RootMpiPNtupleManager::CreateNtuple: wait_buffer() failed."<< G4endl;
    return;
  }
  G4cout << "After wait_buffer with " << fImpi << G4endl;

  tools::uint32 mainNtupleId;
  bool rowWise;
  bool byteSwap;
  unsigned int compression;
  tools::wroot::seek seekDirectory;
  tools::uint32 basketSize;
  bool rowMode;
  std::vector<tools::uint32> basketSizes;
  unsigned int basketEntries;

  if ( ! fImpi->unpack(mainNtupleId) ) {
    G4cerr << "bunpack(byteSwap) failed."<< G4endl;
    return;
  }
  // G4cout << "unpack 1" << G4endl;
  
  if ( ! fImpi->bunpack(rowWise) ) {
    G4cerr << "bunpack(rowWise) failed."<< G4endl;
    return;
  }
  // G4cout << "unpack 1/2" << G4endl;

  if ( ! fImpi->bunpack(byteSwap) ) {
    G4cerr << "bunpack(byteSwap) failed."<< G4endl;
    return;
  }
  // G4cout << "unpack 2" << G4endl;
  
  if ( ! fImpi->unpack(compression) ) {
    G4cerr << "unpack(compression) failed."<< G4endl;
    return;
  }
  // G4cout << "unpack 3" << G4endl;

  if ( ! fImpi->unpack(seekDirectory) ) {
    G4cerr << "unpack(seek) failed."<< G4endl;
    return;
  }
  // G4cout << "unpack 4" << G4endl;
  
  if (rowWise) {
    if ( ! fImpi->unpack(basketSize) ) {
      G4cerr << "unpack(basketSize) failed."<< G4endl;
      return;
    }
  } else {
    if ( ! fImpi->bunpack(rowMode) ) {
      G4cerr << "bpack(rowMode) failed." << G4endl;
      return;
    }      
    if ( ! fImpi->vunpack(basketSizes) ) {
      G4cerr << "vunpack(basketSizes) failed."<< G4endl;
      return;
    }
    if( ! fImpi->unpack(basketEntries) ) {
      G4cerr << "unpack(basketEntries) failed." << G4endl;
      return;
    }
  }
  // G4cout << "unpack 5" << G4endl;

  // G4cout << "rank = " << fMpiRank
  //           << ", verbose = " << verbose
  //           << ", mainNtupleId = " << mainNtupleId 
  //           << ", byteSwap = " << byteSwap
  //           << ", compression = " << compression
  //           << ", seekDirectory = " << seekDirectory
  //           << ", basketSizes.size() = " << basketSizes.size()
  //           << G4endl;

  // Create MPI pntuple
  if ( rowWise ) {
    tools::wroot::mpi_ntuple_row_wise* ntuple
      = new tools::wroot::mpi_ntuple_row_wise(
              mainNtupleId, G4cout, byteSwap, compression, seekDirectory,
              basketSize, ntupleDescription->fNtupleBooking, verbose);
    ntupleDescription->fNtuple = ntuple;
    ntupleDescription->fBasePNtuple = ntuple; 
  } else {
    tools::wroot::mpi_ntuple_column_wise* ntuple
      = new tools::wroot::mpi_ntuple_column_wise(
              mainNtupleId, G4cout, byteSwap, compression, seekDirectory,
              basketSizes, ntupleDescription->fNtupleBooking, 
              rowMode, basketEntries, verbose);
    ntupleDescription->fNtuple = ntuple;
    ntupleDescription->fBasePNtuple = ntuple; 
  }

  ntupleDescription->fIsNtupleOwner = true; 
  ntupleDescription->fImpi = fImpi; 
         // should be not needed 
         // pntuple object is not deleted automatically
  fNtupleVector.push_back(ntupleDescription->fNtuple);  

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()
        ->Message("create from booking", "mpi pntuple", ntupleDescription->fNtupleBooking.name());
#endif
}

//_____________________________________________________________________________
void G4RootMpiPNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuple from ntuple_booking & the buffer from main rank

  // G4cout << "Start G4RootMpiPNtupleManager::CreateNtuplesFromBooking" << G4endl;

  // Do not create ntuples if NtupleVector is not empty
  if ( fNtupleVector.size() ) return;

  for ( auto ntupleDescription : fNtupleDescriptionVector ) {

    // Do not create ntuple if it is inactivated 
    if ( fState.GetIsActivation() && ( ! ntupleDescription->fActivation ) ) continue;
    
    // Do not create ntuple if it already exists
    if ( ntupleDescription->fNtuple ) continue;
    
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from booking", "mpi pntuple", 
                   ntupleDescription->fNtupleBooking.name());
#endif

    // create ntuple
    CreateNtuple(ntupleDescription);

    // finish created ntuple 
    // (nothing to be done ?)
    // FinishTNtuple(ntupleDescription);

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()
        ->Message("create from booking", "mpi pntuple", 
                  ntupleDescription->fNtupleBooking.name());
#endif
  }

}   

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::CreateNtuple(
  const G4String& name, const G4String& title)
{
// Create pntuple description with ntuple_booking

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "pntuple booking", name);
#endif

  // Create ntuple description
  auto index = fNtupleDescriptionVector.size();
  auto ntupleDescription = new G4RootMpiPNtupleDescription();
  ntupleDescription->fMainNtupleRank = fDestinationRank;
        // make this flexible
  fNtupleDescriptionVector.push_back(ntupleDescription);  

  // Save name & title in ntuple booking
  ntupleDescription->fNtupleBooking.set_name(name);
  ntupleDescription->fNtupleBooking.set_title(title);

  fLockFirstId = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << index + fFirstId;
    fState.GetVerboseL2()->Message("create", "pntuple booking", description);
  } 
#endif

  return index + fFirstId;
}                                         

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::CreateNtupleIColumn(
  G4int ntupleId, const G4String& name, std::vector<int>* vector)
{
  return CreateNtupleTColumn<int>(ntupleId, name, vector);
}                                         

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::CreateNtupleFColumn(
  G4int ntupleId, const G4String& name, std::vector<float>* vector)
{
  return CreateNtupleTColumn<float>(ntupleId, name, vector);
}                                         

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::CreateNtupleDColumn(
  G4int ntupleId, const G4String& name, std::vector<double>* vector)
{
  return CreateNtupleTColumn<double>(ntupleId, name, vector);
}                                         

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::CreateNtupleSColumn(
  G4int ntupleId, const G4String& name)
{
  return CreateNtupleTColumn<std::string>(ntupleId, name, nullptr);
}  

//_____________________________________________________________________________
void G4RootMpiPNtupleManager::FinishNtuple(G4int /*ntupleId*/)
{}

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::FillNtupleIColumn(
  G4int ntupleId, G4int columnId, G4int value)
{
  return FillNtupleTColumn<int>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::FillNtupleFColumn(
  G4int ntupleId, G4int columnId, G4float value)
{
  return FillNtupleTColumn<float>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::FillNtupleDColumn(
  G4int ntupleId, G4int columnId, G4double value)
{
  return FillNtupleTColumn<double>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::FillNtupleSColumn(
  G4int ntupleId, G4int columnId, const G4String& value)
{
  return FillNtupleTColumn<std::string>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::AddNtupleRow(G4int ntupleId)
{ 
  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    //G4cout << "Skipping AddNtupleRow for " << ntupleId << G4endl; 
    return false; 
  }  

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("add", "pntuple row", description);
  }  
#endif

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "AddNtupleRow");
  if ( ! ntupleDescription ) return false;
  
  // if(!_ntuple->add_row(_impi,rank_src,tag)) {
  auto result 
    = ntupleDescription->fNtuple
      ->add_row(*fImpi, ntupleDescription->fMainNtupleRank, kTAG_NTUPLE );

  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << "adding row has failed.";
    G4Exception("G4RootMpiPNtupleManager::AddNtupleRow()",
                "Analysis_W002", JustWarning, description);
  }         

#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL3()->Message("add", "pntuple row", description);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::Merge()
{

  // G4cout << "Start G4RootMpiPNtupleManager::Merge" << G4endl;

  auto finalResult = true;

  for ( auto ntupleDescription : fNtupleDescriptionVector) {

    // skip inactivated ntuples
    if ( ! ntupleDescription->fActivation ) continue;
  
    // skip if ntuple was already merged and deleted
    // (this happend when the main rank re-opens a new file for merged data)
    if ( ! ntupleDescription->fNtuple ) continue;
  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()
        ->Message("end_fill", "pntuple", ntupleDescription->fNtupleBooking.name());
    }  
#endif
  
    // if(!_ntuple->end_fill(_impi,rank_src,tag)) {} //important.
    G4cout << "call end_fill " << fImpi 
           << " to rank " << ntupleDescription->fMainNtupleRank 
           << " pntuple: " << ntupleDescription->fNtuple << G4endl;
    auto result 
      = ntupleDescription->fNtuple
        ->end_fill(*fImpi, ntupleDescription->fMainNtupleRank, kTAG_NTUPLE);

    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << " ntuple " << ntupleDescription->fNtupleBooking.name()
                  << "end fill has failed.";
      G4Exception("G4RootMpiPNtupleManager::Merge()",
                  "Analysis_W002", JustWarning, description);
    }

    delete ntupleDescription->fNtuple;
    ntupleDescription->fNtuple = nullptr;
  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) {
      fState.GetVerboseL3()
        ->Message("end_fill", "pntuple", ntupleDescription->fNtupleBooking.name());
    }  
#endif
  }

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::Reset(G4bool deleteNtuple)
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    if ( deleteNtuple ) {
      delete ntupleDescription->fNtuple;
    }  
    ntupleDescription->fNtuple = nullptr;
  }

  fNtupleVector.clear(); 

  return true;
}  

//_____________________________________________________________________________

void  G4RootMpiPNtupleManager::SetActivation(
  G4bool activation)
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    ntupleDescription->fActivation = activation;
  } 
}

//_____________________________________________________________________________

void  G4RootMpiPNtupleManager::SetActivation(
  G4int ntupleId, G4bool activation)
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetActivation");
  if ( ! ntupleDescription ) return;

  ntupleDescription->fActivation = activation;
}

//_____________________________________________________________________________
G4bool  G4RootMpiPNtupleManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "GetActivation");
  if ( ! ntupleDescription ) return false;

  return ntupleDescription->fActivation;
}

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::GetNofNtuples() const
{
  return fNtupleVector.size();
}

//_____________________________________________________________________________
G4int G4RootMpiPNtupleManager::GetNofNtupleBookings() const
{
  return fNtupleDescriptionVector.size();
}

//_____________________________________________________________________________
G4bool G4RootMpiPNtupleManager::IsEmpty() const
{
  return ! fNtupleDescriptionVector.size();
}  

//_____________________________________________________________________________
unsigned int G4RootMpiPNtupleManager::GetBasketSize() const
{ 
  if ( ! fFileManager ) {
    G4String inFunction = "G4RootMpiPNtupleManager::::GetBasketSize";
    G4ExceptionDescription description;
    description << "      " << "File manager must be defined first.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return 0;         
  }

  return fFileManager->GetBasketSize(); 
}
