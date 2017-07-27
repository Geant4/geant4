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
// $Id: G4RootNtupleManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4RootFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"

using namespace G4Analysis;

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state,
                                         G4int nofMainManagers)
 : G4TNtupleManager<tools::wroot::ntuple>(state),
   fCreateMode(G4NtupleCreateMode::kUndefined),
   fFileManager(nullptr),
   fNtupleDirectory(nullptr),
   fMainNtupleManagers()
{
  for ( G4int i=0; i<nofMainManagers; ++i) {
    fMainNtupleManagers.push_back(
      new G4RootMainNtupleManager(this, fState));
  }
}

//_____________________________________________________________________________
G4RootNtupleManager::~G4RootNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::SetCreateMode()
{
// Set create mode if not yet defined

#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("set", "ntuple create mode", "");
#endif

  G4String createMode;
  if ( fCreateMode == G4NtupleCreateMode::kUndefined ) {
    if ( fMainNtupleManagers.size() ) {
      if ( fFileManager->GetNtupleFile(0) ) {
        fCreateMode = G4NtupleCreateMode::kMainAfterOpen;
        createMode = "G4NtupleCreateMode::kMainAfterOpen";
      } else {
        fCreateMode = G4NtupleCreateMode::kMainBeforeOpen;
        createMode = "G4NtupleCreateMode::kMainBeforeOpen";
      }
    }
    else {
      if ( fNtupleDirectory ) {
        fCreateMode = G4NtupleCreateMode::kNoMergeAfterOpen;
        createMode = "G4NtupleCreateMode::kNoMergeAfterOpen";
      } else {
        fCreateMode = G4NtupleCreateMode::kNoMergeBeforeOpen;
        createMode = "G4NtupleCreateMode::kNoMergeBeforeOpen";
      }
    }
  }

#ifdef G4VERBOSE
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("set", "ntuple create mode", createMode);
#endif
}

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtuple(
  G4TNtupleDescription<tools::wroot::ntuple>* ntupleDescription,
  const G4String& name, const G4String& title)
{
  // Set create mode if not yet defined
  SetCreateMode();

  if ( fCreateMode == G4NtupleCreateMode::kNoMergeAfterOpen ) {
    if ( ! fNtupleDirectory ) {
      G4String inFunction = "G4RootNtupleManager::::CreateTNtuple";
      G4ExceptionDescription description;
      description << "      " 
        << "Cannot create ntuple. Ntuple directory does not exist." << G4endl;
      G4Exception(inFunction, "Analysis_W002", JustWarning, description);
      return;
    }

    ntupleDescription->fNtuple
      = new tools::wroot::ntuple(*fNtupleDirectory, name, title);
    ntupleDescription->fIsNtupleOwner = false;  
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);       
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtupleFromBooking(
  G4TNtupleDescription<tools::wroot::ntuple>* ntupleDescription)
{
  if ( fCreateMode == G4NtupleCreateMode::kNoMergeBeforeOpen ) {

    if ( ! fNtupleDirectory ) {
      G4String inFunction = "G4RootNtupleManager::::CreateTNtuple";
      G4ExceptionDescription description;
      description << "      " 
        << "Cannot create ntuple. Ntuple directory does not exist." << G4endl;
      G4Exception(inFunction, "Analysis_W002", JustWarning, description);
      return;
    }
    
    ntupleDescription->fNtuple
      = new tools::wroot::ntuple(
              *fNtupleDirectory, ntupleDescription->fNtupleBooking);

    auto basketSize = fFileManager->GetBasketSize();
    ntupleDescription->fNtuple->set_basket_size(basketSize);
 
    ntupleDescription->fIsNtupleOwner = false;  
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);
  }

  if ( fCreateMode == G4NtupleCreateMode::kMainBeforeOpen ) {
    auto counter = 0;
    for ( auto manager : fMainNtupleManagers ) {
      if ( ! manager->GetNtupleVector().size() ) {
        // Create only once !!
        manager->SetNtupleFile(fFileManager->GetNtupleFile(counter));
        manager->SetNtupleDirectory(fFileManager->GetMainNtupleDirectory(counter++));
        manager->CreateNtuplesFromBooking();
      }
    }
  }
}

//_____________________________________________________________________________
void G4RootNtupleManager::FinishTNtuple(
  G4TNtupleDescription<tools::wroot::ntuple>* ntupleDescription)
{
// Create main ntuples

  if ( fCreateMode == G4NtupleCreateMode::kMainAfterOpen ) {
    auto counter = 0;
    for ( auto manager : fMainNtupleManagers ) {
      auto warn = true;
      manager->SetNtupleFile(fFileManager->GetNtupleFile(counter));
      manager->SetNtupleDirectory(fFileManager->GetMainNtupleDirectory(counter++));
      manager->CreateNtuple(ntupleDescription->fNtupleBooking, warn);
    }
  }
}

//_____________________________________________________________________________
G4bool G4RootNtupleManager::Reset(G4bool deleteNtuple)
{
  G4TNtupleManager<tools::wroot::ntuple> ::Reset(deleteNtuple);
    // this will clear ntuple vector

  if ( fCreateMode == G4NtupleCreateMode::kNoMergeAfterOpen ) {
    // clear also ntuple description vector
    fNtupleDescriptionVector.clear();
  }

  auto finalResult = true;
  for ( auto manager : fMainNtupleManagers ) {
    auto result = manager->Reset(false);
    finalResult = result && finalResult;
  }

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootNtupleManager::Merge()
{
  auto finalResult = true;

  for ( auto manager : fMainNtupleManagers ) {
    auto result = manager->Merge();
    finalResult = result && finalResult;
  }

  return finalResult;
}

//_____________________________________________________________________________
G4RootMainNtupleManager* G4RootNtupleManager::GetMainNtupleManager(G4int index) const
{
  if ( index < 0 || index >= G4int(fMainNtupleManagers.size()) ) {
    G4String inFunction = "G4RootNtupleManager::::GetMainNtupleManager";
    G4ExceptionDescription description;
    description << "      " << "main ntuple manager " << index << " does not exist.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return nullptr;         
  }

  return fMainNtupleManagers[index];
}

//_____________________________________________________________________________
unsigned int G4RootNtupleManager::GetBasketSize() const
{ 
  if ( ! fFileManager ) {
    G4String inFunction = "G4RootNtupleManager::::GetBasketSize";
    G4ExceptionDescription description;
    description << "      " << "File manager must be defined first.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return 0;         
  }

  return fFileManager->GetBasketSize(); 
}

