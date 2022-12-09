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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"

using namespace G4Analysis;
using std::to_string;

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state,
                       const std::shared_ptr<G4NtupleBookingManager>& bookingManger,
                       G4int nofMainManagers, G4int nofFiles,
                       G4bool rowWise, G4bool rowMode)
 : G4TNtupleManager<tools::wroot::ntuple, G4RootFile>(state),
   fRowWise(rowWise),
   fRowMode(rowMode)
{
  for ( G4int i=0; i<nofMainManagers; ++i) {
    auto fileNumber = i;
    if ( (i == 0) && (nofFiles == 0) ) {
      // the main ntuple file will be merged in the default file
      fileNumber = -1;
    }
    fMainNtupleManagers.push_back(
      std::make_shared<G4RootMainNtupleManager>(
        this, bookingManger, rowWise, fileNumber, fState));
  }
}

//
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtupleFromBooking(
  RootNtupleDescription* ntupleDescription)
{
  if (fMainNtupleManagers.size() == 0u) {
    // No merging
    if (ntupleDescription->GetNtuple() != nullptr) {
      Warn("Cannot create ntuple. Ntuple already exists.",
        fkClass, "CreateTNtupleFromBooking");
      return;
    }

    // Create ntuple file from ntuple description
    auto ntupleFile = fFileManager->CreateNtupleFile(ntupleDescription);
    if ( ! ntupleFile ) {
      Warn("Cannot create ntuple. Ntuple file does not exist.",
        fkClass, "CreateTNtupleFromBooking");
      return;
    }

    auto directory = std::get<2>(*ntupleFile);
    ntupleDescription->SetNtuple(
      new tools::wroot::ntuple(
            *directory, ntupleDescription->GetNtupleBooking(), fRowWise));

    auto basketSize = fFileManager->GetBasketSize();
    ntupleDescription->GetNtuple()->set_basket_size(basketSize);

    ntupleDescription->SetIsNtupleOwner(false);
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->GetNtuple());
  }
  else {
    // Merging activated
    for ( const auto& manager : fMainNtupleManagers ) {
      manager->CreateNtuple(ntupleDescription);
    }
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::FinishTNtuple(
  RootNtupleDescription* /*ntupleDescription*/, G4bool /*fromBooking*/)
{
  // nothing to be done
}

//_____________________________________________________________________________
G4bool G4RootNtupleManager::Reset()
{
  G4TNtupleManager<tools::wroot::ntuple, G4RootFile> ::Reset();
    // this will clear ntuple vector

  auto result = true;

  for ( const auto& manager : fMainNtupleManagers ) {
    result &= manager->Reset();
  }

  return result;
}

//_____________________________________________________________________________
void G4RootNtupleManager::Clear()
{
  G4TNtupleManager<tools::wroot::ntuple, G4RootFile> ::Clear();
    // this will clear ntuple vector

  for ( const auto& manager : fMainNtupleManagers ) {
    manager->ClearData();
  }
}

//_____________________________________________________________________________
G4bool G4RootNtupleManager::Merge()
{
  auto result = true;

  for ( const auto& manager : fMainNtupleManagers ) {
    result &= manager->Merge();
  }

  return result;
}

//_____________________________________________________________________________
void  G4RootNtupleManager::SetFileManager(
  const std::shared_ptr<G4RootFileManager>& fileManager)
{
  fFileManager = fileManager;

  for ( const auto& manager : fMainNtupleManagers ) {
    manager->SetFileManager(fileManager);
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::SetNtupleRowWise(G4bool rowWise, G4bool rowMode)
{
// Set rowWise mode and propagate it to main ntuple managers

  fRowWise = rowWise;
  fRowMode = rowMode;

  for (auto& mainNtupleManager : fMainNtupleManagers ) {
    mainNtupleManager->SetRowWise(rowWise);
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::SetNewCycle(G4bool value)
{
  G4TNtupleManager<tools::wroot::ntuple, G4RootFile>::SetNewCycle(value);

  for ( const auto& manager : fMainNtupleManagers ) {
    manager->SetNewCycle(value);
  }
}

//_____________________________________________________________________________
std::shared_ptr<G4RootMainNtupleManager>
G4RootNtupleManager::GetMainNtupleManager(G4int index) const
{
  if ( index < 0 || index >= G4int(fMainNtupleManagers.size()) ) {
    Warn("main ntuple manager " + to_string(index) + " does not exist.",
      fkClass, "GetMainNtupleManager");
    return nullptr;
  }

  return fMainNtupleManagers[index];
}

//_____________________________________________________________________________
unsigned int G4RootNtupleManager::GetBasketSize() const
{
  if ( ! fFileManager ) {
    Warn("File manager must be defined first.", fkClass, "GetBasketSize");
    return 0;
  }

  return fFileManager->GetBasketSize();
}

//_____________________________________________________________________________
unsigned int G4RootNtupleManager::GetBasketEntries() const
{
  if ( ! fFileManager ) {
    Warn("File manager must be defined first.", fkClass, "GetBasketEntries");
    return 0;
  }

  return fFileManager->GetBasketEntries();
}
