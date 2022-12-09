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
 : G4TRNtupleManager<tools::rroot::ntuple>(state)
{}

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
  auto rfileTuple = fFileManager->GetRFile(fileName, isPerThread);
  if (rfileTuple == nullptr) {
    if ( ! fFileManager->OpenRFile(fileName, isPerThread) ) return kInvalidId;
    rfileTuple = fFileManager->GetRFile(fileName, isPerThread);
  }
  auto rfile = std::get<0>(*rfileTuple);

  // Get or open a directory if dirName is defined
  tools::rroot::TDirectory* ntupleDirectory = nullptr;
  if ( ntupleDirectory == nullptr ) {
    // Retrieve directory only once as analysis manager supports only
    // 1 ntuple directory par file)
    if ( ! dirName.empty() ) {
      ntupleDirectory = tools::rroot::find_dir(rfile->dir(), dirName);
      if ( ntupleDirectory != nullptr ) {
        std::get<2>(*rfileTuple) = ntupleDirectory;
      }
      else {
        G4Analysis::Warn(
          "Directory " + dirName + " not found in file " + fileName + ".",
          fkClass, "ReadNtupleImpl");
        return kInvalidId;
      }
    }
  }

  // Get key
  tools::rroot::key* key = nullptr;
  if ( ntupleDirectory != nullptr ) {
    key = ntupleDirectory->find_key(ntupleName);
  }
  else {
    key = rfile->dir().find_key(ntupleName);
  }
  if (key == nullptr) {
    Warn("Key " + ntupleName + " for Ntuple not found in file " + fileName +
      ", directory " + dirName, fkClass, "ReadNtupleImpl");
    return kInvalidId;
  }

  unsigned int size;
  char* charBuffer = key->get_object_buffer(*rfile, size);
  if (charBuffer == nullptr) {
    Warn("Cannot get data buffer for Ntuple " + ntupleName +
         " in file " + fileName, fkClass, "ReadNtupleImpl");
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
  auto ntuple = ntupleDescription->fNtuple;
  auto ntupleBinding = ntupleDescription->fNtupleBinding;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {

    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      Warn("Ntuple initialization failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  auto next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      Warn("Ntuple get_row() failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
  }

  return next;
}
