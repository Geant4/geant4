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
// $Id$
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#include "G4RootPNtupleManager.hh"
#include "G4RootMainNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/ntuple"

// mutex in a file scope
namespace {
  //Mutex to lock master manager when adding ntuple row
  G4Mutex addRowMutex = G4MUTEX_INITIALIZER;
  G4Mutex endFillMutex = G4MUTEX_INITIALIZER;
}

//_____________________________________________________________________________
G4RootPNtupleManager::G4RootPNtupleManager(G4RootMainNtupleManager* main,
                                           const G4AnalysisManagerState& state)
 : G4BaseNtupleManager(state),
   fCreateMode(G4PNtupleCreateMode::kUndefined),
   fMainNtupleManager(main),
   fNtupleVector()
{}

//_____________________________________________________________________________
G4RootPNtupleManager::~G4RootPNtupleManager()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    delete ntupleDescription;
  }   
}

//
// private functions
//

//_____________________________________________________________________________
G4RootPNtupleDescription*  
G4RootPNtupleManager::GetNtupleDescriptionInFunction(
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootPNtupleManager::";
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
tools::wroot::base_pntuple*  
 G4RootPNtupleManager::GetNtupleInFunction(
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto ntupleDescription = GetNtupleDescriptionInFunction(id, functionName);
  if ( ! ntupleDescription ) return nullptr;

  if ( ! ntupleDescription->fBasePNtuple ) {
    if ( warn ) {
      G4String inFunction = "G4RootPNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      "  << "ntupleId " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return nullptr;
  }  
  return ntupleDescription->fBasePNtuple;
}

//_____________________________________________________________________________
tools::wroot::ntuple*
G4RootPNtupleManager::GetMainNtupleInFunction(
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto& mainNtupleVector
     = fMainNtupleManager->GetNtupleVector();
  
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(mainNtupleVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4RootPNtupleManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "main ntuple " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return nullptr;         
  }
  
  return mainNtupleVector[index];
}


//
// protected functions
//

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtuple(G4RootPNtupleDescription* ntupleDescription,
                                        tools::wroot::ntuple* mainNtuple)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()
        ->Message("create from main", "pntuple", mainNtuple->name());
#endif

  auto rfile = fMainNtupleManager->GetNtupleFile();
  // auto directory = fMainNtupleManager->GetNtupleDirectory();

  // Get parameters from main ntuple
  auto mainBranch = mainNtuple->get_row_wise_branch();
  auto rowWise = mainBranch ? true : false;

  ntupleDescription->fFile = rfile.get();
  mainNtuple->get_branches(ntupleDescription->fMainBranches);

  G4bool verbose = true;
  if ( rowWise ) {
    tools::wroot::mt_ntuple_row_wise* mtNtuple 
      = new tools::wroot::mt_ntuple_row_wise(
              G4cout, rfile->byte_swap(), rfile->compression(),
              mainNtuple->dir().seek_directory(),
              *mainBranch, mainBranch->basket_size(), 
              ntupleDescription->fNtupleBooking, verbose);

    ntupleDescription->fNtuple
      = static_cast<tools::wroot::imt_ntuple*>(mtNtuple);
    ntupleDescription->fBasePNtuple
      = static_cast<tools::wroot::base_pntuple*>(mtNtuple);
  } 
  else {
    std::vector<tools::uint32> basketSizes;
    tools_vforcit(tools::wroot::branch*, ntupleDescription->fMainBranches, it) {
      basketSizes.push_back((*it)->basket_size());
    }
   
    tools::wroot::mt_ntuple_column_wise* mtNtuple =
      new tools::wroot::mt_ntuple_column_wise(
            G4cout, rfile->byte_swap(), rfile->compression(),
            mainNtuple->dir().seek_directory(),
            ntupleDescription->fMainBranches, basketSizes,
            ntupleDescription->fNtupleBooking, verbose);
    
    ntupleDescription->fNtuple
      = static_cast<tools::wroot::imt_ntuple*>(mtNtuple);
    ntupleDescription->fBasePNtuple
      = static_cast<tools::wroot::base_pntuple*>(mtNtuple);
  }

  ntupleDescription->fIsNtupleOwner = true;  
  //        // pntuple object is not deleted automatically
  fNtupleVector.push_back(ntupleDescription->fNtuple);  


#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()
        ->Message("create from main", "pntuple", mainNtuple->name());
#endif
}

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtuplesFromMain()
{
// Create ntuple from booking (if not yet done) and main ntuple 
// This function is called from G4AnalysisManager::OpenFile.

  if ( fCreateMode == G4PNtupleCreateMode::kUndefined ) {
    if ( fNtupleDescriptionVector.size() ) {
      fCreateMode = G4PNtupleCreateMode::kSlaveBeforeOpen;
      // G4cout << "Create mode: kSlaveBeforeOpen" <<  G4endl;
    } else {
      fCreateMode = G4PNtupleCreateMode::kSlaveAfterOpen;
      // G4cout << "Create mode: kSlaveAfterOpen" <<  G4endl;
    }
  }

  if ( fCreateMode == G4PNtupleCreateMode::kSlaveAfterOpen ) {
    // ntuples are not yet booked
    // G4cout << "Ntuples are not yet booked ?" <<  G4endl;
    return;
  }

  auto& mainNtupleVector
     = fMainNtupleManager->GetNtupleVector();

  G4int lcounter = 0;
  for ( auto mainNtuple : mainNtupleVector ) {
    
    auto& ntupleDescription 
      = fNtupleDescriptionVector[lcounter++];
    CreateNtuple(ntupleDescription, mainNtuple);
  }
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::CreateNtuple(
  const G4String& name, const G4String& title)
{
// Create pntuple description with ntuple_booking

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "pntuple booking", name);
#endif

  // Set create mode if not yet defined
  if ( fCreateMode == G4PNtupleCreateMode::kUndefined ) {
    if ( fMainNtupleManager->GetNtupleFile() ) {
      fCreateMode = G4PNtupleCreateMode::kSlaveAfterOpen;
    } else {
      fCreateMode = G4PNtupleCreateMode::kSlaveBeforeOpen;
    }
  }

  // Create ntuple description
  auto index = fNtupleDescriptionVector.size();
  auto ntupleDescription = new G4RootPNtupleDescription();
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
G4int G4RootPNtupleManager::CreateNtupleIColumn(
  G4int ntupleId, const G4String& name, std::vector<int>* vector)
{
  return CreateNtupleTColumn<int>(ntupleId, name, vector);
}                                         

//_____________________________________________________________________________
G4int G4RootPNtupleManager::CreateNtupleFColumn(
  G4int ntupleId, const G4String& name, std::vector<float>* vector)
{
  return CreateNtupleTColumn<float>(ntupleId, name, vector);
}                                         


//_____________________________________________________________________________
G4int G4RootPNtupleManager::CreateNtupleDColumn(
  G4int ntupleId, const G4String& name, std::vector<double>* vector)
{
  return CreateNtupleTColumn<double>(ntupleId, name, vector);
}                                         

//_____________________________________________________________________________
G4int G4RootPNtupleManager::CreateNtupleSColumn(
  G4int ntupleId, const G4String& name)
{
  return CreateNtupleTColumn<std::string>(ntupleId, name, nullptr);
}  

//_____________________________________________________________________________
void G4RootPNtupleManager::FinishNtuple(G4int ntupleId)
{ 
// create ntuple if file was open

  // if ( fMainNtupleManager->GetNtupleFile() ) {
  if ( fCreateMode == G4PNtupleCreateMode::kSlaveAfterOpen ) {
    auto ntupleDescription 
      = GetNtupleDescriptionInFunction(ntupleId, "FinishNtuple");
    if ( ! ntupleDescription ) return;

    auto mainNtuple = GetMainNtupleInFunction(ntupleId, "FinishNtuple");
    if ( ! mainNtuple ) return;
    
    CreateNtuple(ntupleDescription, mainNtuple);
  }
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::FillNtupleIColumn(
  G4int ntupleId, G4int columnId, G4int value)
{
  return FillNtupleTColumn<int>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::FillNtupleFColumn(
  G4int ntupleId, G4int columnId, G4float value)
{
  return FillNtupleTColumn<float>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::FillNtupleDColumn(
  G4int ntupleId, G4int columnId, G4double value)
{
  return FillNtupleTColumn<double>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::FillNtupleSColumn(
  G4int ntupleId, G4int columnId, const G4String& value)
{
  return FillNtupleTColumn<std::string>(ntupleId, columnId, value);
}                                         

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::AddNtupleRow(G4int ntupleId)
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
  
  G4AutoLock lock(&addRowMutex);
  mutex toolsLock(lock);
  auto result 
    = ntupleDescription->fNtuple
      ->add_row(toolsLock, *ntupleDescription->fFile);

  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << " ntupleId " << ntupleId 
                << "adding row has failed.";
    G4Exception("G4RootPNtupleManager::AddNtupleRow()",
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
G4bool G4RootPNtupleManager::Merge()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector) {

    // skip inactivated ntuples
    if ( ! ntupleDescription->fActivation ) continue;
  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) {
      fState.GetVerboseL4()
        ->Message("merge", "pntuple", ntupleDescription->fNtupleBooking.name());
    }  
#endif
  
    G4AutoLock lock(&endFillMutex);
    mutex toolsLock(lock);
    auto result 
      = ntupleDescription->fNtuple
        ->end_fill(toolsLock, *ntupleDescription->fFile);

    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << " ntuple " << ntupleDescription->fNtupleBooking.name()
                  << "end fill has failed.";
      G4Exception("G4RootPNtupleManager::Merge()",
                  "Analysis_W002", JustWarning, description);
    }

    delete ntupleDescription->fNtuple;
    ntupleDescription->fNtuple = nullptr;
  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) {
      fState.GetVerboseL3()
        ->Message("merge", "pntuple", ntupleDescription->fNtupleBooking.name());
    }  
#endif
  }
  return true;

}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::Reset(G4bool deleteNtuple)
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

void  G4RootPNtupleManager::SetActivation(
  G4bool activation)
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    ntupleDescription->fActivation = activation;
  } 
}

//_____________________________________________________________________________

void  G4RootPNtupleManager::SetActivation(
  G4int ntupleId, G4bool activation)
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetActivation");
  if ( ! ntupleDescription ) return;

  ntupleDescription->fActivation = activation;
}

//_____________________________________________________________________________
G4bool  G4RootPNtupleManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "GetActivation");
  if ( ! ntupleDescription ) return false;

  return ntupleDescription->fActivation;
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::GetNofNtuples() const
{
  return fNtupleVector.size();
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::GetNofNtupleBookings() const
{
  return fNtupleDescriptionVector.size();
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::IsEmpty() const
{
  return ! fNtupleDescriptionVector.size();
}  

