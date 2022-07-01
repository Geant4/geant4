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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4Hdf5RNtupleManager.hh"
#include "G4Hdf5RFileManager.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4Hdf5RNtupleManager::G4Hdf5RNtupleManager(const G4AnalysisManagerState& state)
 : G4TRNtupleManager<toolx::hdf5::ntuple>(state)
{}

//
// private methods
//

//_____________________________________________________________________________
G4int G4Hdf5RNtupleManager::ReadNtupleImpl(const G4String& ntupleName,
                                           const G4String& fileName,
                                           const G4String& dirName,
                                           G4bool isUserFileName)
{
  Message(kVL4, "read", "ntuple", ntupleName);

  // Ntuples are saved in files per thread
  // but apply thethe thread suffix only if fileName is not provided explicitly
  G4String fullFileName = fileName;
  if ( ! isUserFileName ) {
    fullFileName = fFileManager->GetFullFileName();
  }

  // Get directory
  auto directory = fFileManager->GetNtupleRDirectory(fullFileName, dirName, false);
  if ( directory < 0 ) return kInvalidId;

  // Create ntuple
  auto rntuple = new toolx::hdf5::ntuple(G4cout, directory, ntupleName);
  auto rntupleDescription = new G4TRNtupleDescription<toolx::hdf5::ntuple>(rntuple);
  auto id = SetNtuple(rntupleDescription);

  Message(kVL2, "read", "ntuple", ntupleName, id > kInvalidId);

  return id;
}

//_____________________________________________________________________________
G4bool G4Hdf5RNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<toolx::hdf5::ntuple>* ntupleDescription)
{
  auto ntuple = ntupleDescription->fNtuple;

  auto isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {
    auto ntupleBinding = ntupleDescription->fNtupleBinding;
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      Warn("Ntuple initialization failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
    ntupleDescription->fIsInitialized = true;
  }

  if ( ! ntuple->get_row() ) {
    Warn( "Ntuple get_row() failed !!", fkClass, "GetTNtupleRow");
    return false;
  }

  return true;
}
