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

  // Get ntuple booking
  const auto& ntupleBooking = ntupleDescription->GetNtupleBooking();

  Message(kVL4, "create", "main ntuple", ntupleBooking.name());

  // Create ntuple
  auto ntuple = new tools::wroot::ntuple(*std::get<2>(*ntupleFile), ntupleBooking, fRowWise);
         // ntuple object is deleted automatically when closing a file
  auto basketSize = fNtupleBuilder->GetBasketSize();
  ntuple->set_basket_size(basketSize);

  fNtupleVector.push_back(ntuple);
  fNtupleDescriptionVector.push_back(std::make_pair(ntupleDescription, ntupleFile));

  Message(kVL3, "create", "main ntuple", ntupleBooking.name());
}

//_____________________________________________________________________________
G4bool G4RootMainNtupleManager::Merge()
{
  std::size_t counter = 0;

  for ( auto ntuple : fNtupleVector ) {
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
// Create ntuples from booking.
// This function is triggered from workers at new cycle.

  for (auto [ntupleDescription, ntupleFile] : fNtupleDescriptionVector) {

    // Get ntuple booking
    const auto& ntupleBooking = ntupleDescription->GetNtupleBooking();

    Message(kVL4, "create", "main ntuple", ntupleBooking.name());

    // Create ntuple
    auto ntuple = new tools::wroot::ntuple(*std::get<2>(*ntupleFile), ntupleBooking, fRowWise);
         // ntuple object is deleted automatically when closing a file
    auto basketSize = fNtupleBuilder->GetBasketSize();
    ntuple->set_basket_size(basketSize);

    fNtupleVector.push_back(ntuple);

    Message(kVL3, "create", "main ntuple", ntupleBooking.name());
  }

  SetNewCycle(false);
}
