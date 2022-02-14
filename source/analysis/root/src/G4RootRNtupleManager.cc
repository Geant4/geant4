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

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#include "G4RootRNtupleManager.hh"
#include "G4RootRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/rroot/file"
#include "tools/rroot/rall"
#include "tools/rroot/streamers"
#include "tools/rroot/fac"
#include "tools/rroot/tree"

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootRNtupleManager::G4RootRNtupleManager(const G4AnalysisManagerState& state)
 : G4VRNtupleManager(state),
   fNtupleVector()
{
}

//
// private methods
//

//_____________________________________________________________________________
G4int G4RootRNtupleManager::ReadNtupleImpl(const G4String& ntupleName,
                                           const G4String& fileName,
                                           const G4String& dirName,
                                           G4bool isUserFileName)
{
  Message(kVL4, "read", "ntuple", ntupleName);

  // Ntuples are saved per thread
  // but do not apply the thread suffix if fileName is provided explicitly
  auto isPerThread = true;
  if ( isUserFileName ) isPerThread = false;

  // Get or open a file
  auto rfile = fFileManager->GetRFile(fileName, isPerThread);
  if ( ! rfile ) {
    if ( ! fFileManager->OpenRFile(fileName, isPerThread) ) return kInvalidId;
    rfile = fFileManager->GetRFile(fileName, isPerThread);
  }

  // Get key
  tools::rroot::key* key = nullptr;
  tools::rroot::TDirectory* newDir = nullptr;
  if ( ! dirName.empty() ) {
    newDir = tools::rroot::find_dir(rfile->dir(), dirName);
    if ( newDir ) {
      key = newDir->find_key(ntupleName);
    }
    else {
      G4Analysis::Warn(
        "Directory " + dirName + " not found in file " + fileName + ".",
        fkClass, "ReadNtupleImpl");
      return kInvalidId;
    }
  }
  else {
    key = rfile->dir().find_key(ntupleName);
  }

  if ( ! key ) {
    Warn("Key " + ntupleName + " for Ntuple not found in file " + fileName +
      ", directory " + dirName, fkClass, "ReadNtupleImpl");
    delete newDir;
    return kInvalidId;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(*rfile, size);
  if ( ! charBuffer ) {
    Warn("Cannot get data buffer for Ntuple " + ntupleName +
         " in file " + fileName, fkClass, "ReadNtupleImpl");
    delete newDir;
    return kInvalidId;
  }

  auto verbose = false;
  auto buffer
    = new tools::rroot::buffer(G4cout, rfile->byte_swap(), size, charBuffer,
                               key->key_length(), verbose);
  buffer->set_map_objs(true);

  auto fac = new tools::rroot::fac(G4cout);

  auto tree = new tools::rroot::tree(*rfile, *fac);
  if ( ! tree->stream(*buffer) ) {
    Warn("TTree streaming failed for Ntuple " + ntupleName +
         " in file " + fileName,
         fkClass, "ReadNtupleImpl");

    delete buffer;
    delete tree;
    delete newDir;
    return kInvalidId;
  }

  auto rntuple  = new tools::rroot::ntuple(*tree); //use the flat ntuple API.
  auto rntupleDescription = new G4TRNtupleDescription<tools::rroot::ntuple>(rntuple);

  auto id = SetNtuple(rntupleDescription);

  Message(kVL2, "read", "ntuple", ntupleName, id > kInvalidId);

  return id;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<tools::rroot::ntuple>* ntupleDescription)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple I column", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleIColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple I colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleFColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4float>& vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple F column of vector", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleFColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple F colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::SetNtupleDColumn(G4int ntupleId, 
                                             const G4String& columnName, 
                                             std::vector<G4double>& vector)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL4()->Message("set", "ntuple D column of vector", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "SetNtupleDColumn");
  if ( ! ntupleDescription )  return false;   
  
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
  ntupleBinding->add_column(columnName, vector);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId << " " << columnName;  
    fState.GetVerboseL2()->Message("set", "ntuple D colum", description, true);
  }  
#endif

  return true;
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetNtupleRow()
{ 
  return GetNtupleRow(fFirstId);                                            
}

//_____________________________________________________________________________
G4bool G4RootRNtupleManager::GetNtupleRow(G4int ntupleId)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("get", "ntuple row", description);
  }  
#endif

  G4RootRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtupleRow");
  if ( ! ntupleDescription )  return false;   
  
  tools::rroot::ntuple* ntuple = ntupleDescription->fNtuple;
  tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      G4ExceptionDescription description;
      description << " ntupleId " << ntupleId;  
      fState.GetVerboseL4()->Message("initialize", "ntuple", description);
    }  
#endif
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      Warn("Ntuple initialization failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  G4bool next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      Warn("Ntuple get_row() failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL2()->Message("get", "ntuple row", description, true);
  }  
#endif

  return next;
}
