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

#include "G4Hdf5NtupleManager.hh"
#include "G4Hdf5FileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4UnitsTable.hh"

#include "tools/ntuple_booking"

using namespace G4Analysis;

//_____________________________________________________________________________
G4Hdf5NtupleManager::G4Hdf5NtupleManager(const G4AnalysisManagerState& state)
 : G4TNtupleManager<toolx::hdf5::ntuple, G4Hdf5File>(state)
{}

//
// private methods
//
//_____________________________________________________________________________
void G4Hdf5NtupleManager::CreateTNtuple(
  Hdf5NtupleDescription* ntupleDescription, G4bool warn)
{
// Ntuple will be created from ntupleDescription if file is open,
// return with or without warning otherwise

  // Get ntuple file from ntuple description
  auto ntupleFile = ntupleDescription->GetFile();
  if (! ntupleFile) {
    ntupleFile = fFileManager->GetFile();
  }

  // Check file
  if ( ! ntupleFile ) {
    if (warn) {
      Warn( "Cannot create ntuple. Ntuple file does not exist.",
        fkClass, "CreateTNtuple");
    }
    return;
  }

  auto directory = std::get<2>(*ntupleFile);
  auto basketSize = fFileManager->GetBasketSize();
  // auto compressionLevel = fState.GetCompressionLevel();
  auto compressionLevel = 0;

  // Update ntuple name if cycle >0
  auto ntupleBooking = ntupleDescription->GetNtupleBooking();
  if (GetCycle() > 0) {
    auto newNtupleName = ntupleDescription->GetNtupleBooking().name();
    newNtupleName.append("_v");
    newNtupleName.append(std::to_string(GetCycle()));
    ntupleBooking.set_name(newNtupleName);
  }

  // create ntuple
  ntupleDescription->SetNtuple(
    new toolx::hdf5::ntuple(
          G4cout, directory, ntupleBooking, compressionLevel, basketSize));

  fNtupleVector.push_back(ntupleDescription->GetNtuple());
}

//_____________________________________________________________________________
void G4Hdf5NtupleManager::CreateTNtupleFromBooking(
  Hdf5NtupleDescription* ntupleDescription)
{
  // Create file if file name per object is set
  if (ntupleDescription->GetFileName().size() != 0u) {
    fFileManager->CreateNtupleFile(ntupleDescription);
  }

  // Create ntuple from booking and print a warning if file is not open
  CreateTNtuple(ntupleDescription, true);
}

//_____________________________________________________________________________
void G4Hdf5NtupleManager::FinishTNtuple(
  Hdf5NtupleDescription* ntupleDescription,
  G4bool /*fromBooking*/)
{
  if (ntupleDescription->GetNtuple() == nullptr) {
    // Create ntuple from booking if file is open, do nothing otherwise
    CreateTNtuple(ntupleDescription, false);
  }

  fFileManager->LockDirectoryNames();
}

