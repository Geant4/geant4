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

//_____________________________________________________________________________
G4RootNtupleManager::G4RootNtupleManager(const G4AnalysisManagerState& state,
                       std::shared_ptr<G4NtupleBookingManager> bookingManger,
                       G4int nofMainManagers, G4int nofFiles,
                       G4bool rowWise, G4bool rowMode)
 : G4TNtupleManager<tools::wroot::ntuple, G4RootFile>(state),
   fFileManager(nullptr),
   fMainNtupleManagers(),
   fNtupleFile(nullptr),
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
        this, bookingManger.get(), rowWise, fileNumber, fState));
  }
}

//_____________________________________________________________________________
G4RootNtupleManager::~G4RootNtupleManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
void G4RootNtupleManager::CreateTNtupleFromBooking(
  RootNtupleDescription* ntupleDescription)
{
  if ( ! fMainNtupleManagers.size() ) {
    // No merging
    if ( ntupleDescription->fNtuple ) {
      G4String inFunction = "G4RootNtupleManager::::CreateTNtupleFromBooking";
      G4ExceptionDescription description;
      description << "Cannot create ntuple. Ntuple already exists." << G4endl;
      G4Exception(inFunction, "Analysis_W002", JustWarning, description);
      return;
    }        
  
    // Create ntuple file from ntuple description
    auto ntupleFile = fFileManager->CreateNtupleFile(ntupleDescription);
    if ( ! ntupleFile ) {
      G4String inFunction = "G4RootNtupleManager::::CreateTNtupleFromBooking";
      G4ExceptionDescription description;
      description << "Cannot create ntuple. Ntuple file does not exist." << G4endl;
      G4Exception(inFunction, "Analysis_W002", JustWarning, description);
      return;
    }
    
    auto directory = std::get<2>(*ntupleFile);
    ntupleDescription->fNtuple
      = new tools::wroot::ntuple(
              *directory, ntupleDescription->fNtupleBooking, fRowWise);
    
    auto basketSize = fFileManager->GetBasketSize();
    ntupleDescription->fNtuple->set_basket_size(basketSize);
    
    ntupleDescription->fIsNtupleOwner = false;  
           // ntuple object is deleted automatically when closing a file
    fNtupleVector.push_back(ntupleDescription->fNtuple);
  } 
  else {
    // Merging activated  
    for ( auto manager : fMainNtupleManagers ) {
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
G4bool G4RootNtupleManager::Reset(G4bool deleteNtuple)
{
  G4TNtupleManager<tools::wroot::ntuple, G4RootFile> ::Reset(deleteNtuple);
    // this will clear ntuple vector

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
void  G4RootNtupleManager::SetFileManager(std::shared_ptr<G4RootFileManager> fileManager)
{ 
  fFileManager = fileManager; 

  for ( auto mainNtupleManager : fMainNtupleManagers) {
    mainNtupleManager->SetFileManager(fileManager);
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
std::shared_ptr<G4RootMainNtupleManager>
G4RootNtupleManager::GetMainNtupleManager(G4int index) const
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

//_____________________________________________________________________________
unsigned int G4RootNtupleManager::GetBasketEntries() const
{ 
  if ( ! fFileManager ) {
    G4String inFunction = "G4RootNtupleManager::::GetBasketEntries";
    G4ExceptionDescription description;
    description << "      " << "File manager must be defined first.";
    G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    return 0;         
  }

  return fFileManager->GetBasketEntries(); 
}
