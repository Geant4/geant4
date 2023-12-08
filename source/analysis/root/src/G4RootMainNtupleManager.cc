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
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#include "G4RootMainNtupleManager.hh"
#include "G4RootFileManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/ntuple"

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootMainNtupleManager::G4RootMainNtupleManager(
                G4RootNtupleManager* ntupleBuilder,
                std::shared_ptr<G4NtupleBookingManager> bookingManager,
                G4bool rowWise,
                G4int fileNumber,
                const G4AnalysisManagerState& state)
 : G4BaseAnalysisManager(state),
   fNtupleBuilder(ntupleBuilder),
   fBookingManager(std::move(bookingManager)),
   fRowWise(rowWise),
   fFileNumber(fileNumber)
{}

//
// private functions
//

//_____________________________________________________________________________
G4int G4RootMainNtupleManager::CreateNtupleFromBooking(
  G4NtupleBooking* g4NtupleBooking, std::shared_ptr<G4RootFile> ntupleFile)
{
  // Do not create ntuple if it was deleted
  if ( g4NtupleBooking->GetDeleted()) {
    // G4cout << "Ntuple " << g4NtupleBooking->fNtupleId << " was deleted" << G4endl;
    return G4Analysis::kInvalidId;
  }

  // The ntuple index
  auto index = g4NtupleBooking->fNtupleId - fFirstId;

  // Do not create ntuple if it already exists
  if ( (index < G4int(fNtupleVector.size())) && (fNtupleVector[index] != nullptr) ) {
    // G4cout << "Ntuple " << g4NtupleBooking->fNtupleId << " already exists" << G4endl;
    return G4Analysis::kInvalidId;
  }

  // Get ntuple booking
  const auto& ntupleBooking = g4NtupleBooking->fNtupleBooking;

  Message(kVL4, "create", "main ntuple", ntupleBooking.name());

  // Allocate the ntuple vector element(s) if needed
  while ( index >= G4int(fNtupleVector.size()) ) {
    fNtupleVector.push_back(nullptr);
  }

  // Create ntuple
  auto ntuple = new tools::wroot::ntuple(*std::get<2>(*ntupleFile), ntupleBooking, fRowWise);
         // ntuple object is deleted automatically when closing a file
  auto basketSize = fNtupleBuilder->GetBasketSize();
  ntuple->set_basket_size(basketSize);

  // Save ntuple & ntuple description pair in vectors
  fNtupleVector[index] = ntuple;

  Message(kVL3, "create", "main ntuple", ntupleBooking.name());

  return index;
}

//
// protected functions
//

//_____________________________________________________________________________
void G4RootMainNtupleManager::CreateNtuple(RootNtupleDescription* ntupleDescription,
                                           G4bool warn)
{
// Create ntuple from booking if file was open

  // Get/Create main ntuple file
  auto ntupleFile = fFileManager->CreateNtupleFile(ntupleDescription, fFileNumber);
  if ( ntupleFile == nullptr ) {
    if ( warn ) {
      Warn("Ntuple file must be defined first.\n"
           "Cannot create main ntuple.",
           fkClass, "CreateNtuple");
    }
    return;
  }

  // Create ntuple from g4 booking
  auto g4NtupleBooking = ntupleDescription->GetG4NtupleBooking();
  auto index = CreateNtupleFromBooking(g4NtupleBooking, ntupleFile);

  // Return if ntuple was not created
  if (index == G4Analysis::kInvalidId) return;

  // Allocate the ntuple description pair vector element(s) if needed
  while ( index >= G4int(fNtupleDescriptionVector.size()) ) {
    fNtupleDescriptionVector.push_back(std::make_pair(nullptr, nullptr));
  }

  // Save ntuple description pair in vectors
  fNtupleDescriptionVector[index] = std::make_pair(ntupleDescription, ntupleFile);
}

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Delete(G4int id)
{
  if (fNtupleVector.empty()) {
    // Main ntuples are delete with each new cycle
    return true;
  }

  // Proceed with deleting if vector is not empty

  Message(kVL4, "delete", "main ntuple ntupleId: " + to_string(id));

  // Get ntuple description
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleVector.size()) ) {
    G4Analysis::Warn("Main ntuple " + to_string(id) + " does not exist.",
      fkClass, "Delete");
    return false;
  }

  // Delete main ntuple and update ntuple vector
  delete fNtupleVector[index];
  fNtupleVector[index] = nullptr;

  Message(kVL3, "delete", "main ntuple ntupleId: " + to_string(id));

  return true;
}

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Merge()
{
  std::size_t counter = 0;

  for ( auto ntuple : fNtupleVector ) {
    // skip deleted ntuples
    if (ntuple == nullptr) continue;

    ntuple->merge_number_of_entries();

    // Notify ntuple description that file is not empty
    if (ntuple->entries() != 0u) {
      auto ntupleDescription = fNtupleDescriptionVector.at(counter).first;
      ntupleDescription->SetHasFill(true);
    }
    ++counter;
  }

  return true;
}

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Reset()
{
  // The ntuples will be recreated with new cycle or new open file.
  // Ntuple objects are deleted automatically when closing a file

  fNtupleVector.clear();

  return true;
}

//_____________________________________________________________________________
void G4RootMainNtupleManager::ClearData()
{
  fNtupleDescriptionVector.clear();
  fNtupleVector.clear();

  Message(G4Analysis::kVL2, "clear", "main ntuples");
}

//_____________________________________________________________________________
std::shared_ptr<G4RootFile>
G4RootMainNtupleManager::GetNtupleFile(RootNtupleDescription* ntupleDescription) const
{
  auto perThread = false;
  return fFileManager->GetNtupleFile(ntupleDescription, perThread, fFileNumber);
}

//
// public functions
//

//_____________________________________________________________________________
void G4RootMainNtupleManager::CreateNtuplesFromBooking()
{
// Create ntuples from booking (without creating ntuple description)
// This function is triggered from workers at new cycle.

  for (auto [ntupleDescription, ntupleFile] : fNtupleDescriptionVector) {
    CreateNtupleFromBooking(ntupleDescription->GetG4NtupleBooking(), ntupleFile);
  }

  SetNewCycle(false);
}
