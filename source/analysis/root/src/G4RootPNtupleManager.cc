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

#include "G4NtupleBookingManager.hh"
#include "G4RootPNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/ntuple"

using namespace G4Analysis;
using std::to_string;

// mutex in a file scope
namespace {

//Mutex to lock master manager when adding ntuple row and ending fill
G4Mutex pntupleMutex = G4MUTEX_INITIALIZER;
//Mutex to lock master manager when createing main ntuples at new cycle
G4Mutex createMainMutex = G4MUTEX_INITIALIZER;

//_____________________________________________________________________________
void NotExistWarning(const G4String& what, G4int id,
                     std::string_view className,
                     std::string_view functionName)
{
  Warn(what + " id= " + to_string(id) + " does not exist.",
    className, functionName);

}

}

//_____________________________________________________________________________
G4RootPNtupleManager::G4RootPNtupleManager(const G4AnalysisManagerState& state,
                        std::shared_ptr<G4NtupleBookingManager> bookingManger,
                        std::shared_ptr<G4RootMainNtupleManager> main,
                        G4bool rowWise, G4bool rowMode)
 : G4BaseNtupleManager(state),
   fBookingManager(std::move(bookingManger)),
   fMainNtupleManager(std::move(main)),
   fRowWise(rowWise),
   fRowMode(rowMode)
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
  G4int id, std::string_view functionName, G4bool warn) const
{
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleDescriptionVector.size()) ) {
    if ( warn) {
      NotExistWarning("ntuple description", id, fkClass, functionName);
    }
    return nullptr;
  }

  return fNtupleDescriptionVector[index];
}

//_____________________________________________________________________________
tools::wroot::base_pntuple*
 G4RootPNtupleManager::GetNtupleInFunction(
  G4int id, std::string_view functionName, G4bool warn) const
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(id, functionName);
  if (ntupleDescription == nullptr) return nullptr;

  if (ntupleDescription->GetBasePNtuple() == nullptr) {
    if ( warn ) {
      NotExistWarning("ntuple", id, fkClass, functionName);
    }
    return nullptr;
  }
  return ntupleDescription->GetBasePNtuple();
}

//_____________________________________________________________________________
tools::wroot::ntuple*
G4RootPNtupleManager::GetMainNtupleInFunction(
  G4int id, std::string_view functionName, G4bool warn) const
{
  auto& mainNtupleVector = fMainNtupleManager->GetNtupleVector();

  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(mainNtupleVector.size()) ) {
    if ( warn) {
      NotExistWarning("main ntuple", id, fkClass, functionName);
    }
    return nullptr;
  }

  return mainNtupleVector[index];
}

//
// protected functions
//

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtupleFromMain(
                             G4RootPNtupleDescription* ntupleDescription,
                             tools::wroot::ntuple* mainNtuple)
{
  Message(kVL4, "create from main", "pntuple", mainNtuple->name());

  auto file = fMainNtupleManager->GetNtupleFile(&ntupleDescription->GetDescription());
  if ( ! file ) {
    Warn("Cannot create pntuple. Main ntuple file does not exist.",
      fkClass, "CreateNtupleFromMain");
    return;
  }

  ntupleDescription->GetDescription().SetFile(file);

  // Get parameters from ntupleDescription
  mainNtuple->get_branches(ntupleDescription->GetMainBranches());

  auto rfile = std::get<0>(*file);
  G4bool verbose = true;
  if ( fRowWise ) {
    auto mainBranch = mainNtuple->get_row_wise_branch();
    auto mtNtuple
      = new tools::wroot::mt_ntuple_row_wise(
              G4cout, rfile->byte_swap(), rfile->compression(),
              mainNtuple->dir().seek_directory(),
              *mainBranch, mainBranch->basket_size(),
              ntupleDescription->GetDescription().GetNtupleBooking(), verbose);

    ntupleDescription->SetNtuple(
      static_cast<tools::wroot::imt_ntuple*>(mtNtuple));
    ntupleDescription->SetBasePNtuple(
      static_cast<tools::wroot::base_pntuple*>(mtNtuple));
  }
  else {
    std::vector<tools::uint32> basketSizes;
    tools_vforcit(tools::wroot::branch*, ntupleDescription->GetMainBranches(), it) {
      basketSizes.push_back((*it)->basket_size());
    }
    auto basketEntries = fMainNtupleManager->GetBasketEntries();

    auto mtNtuple =
      new tools::wroot::mt_ntuple_column_wise(
            G4cout, rfile->byte_swap(), rfile->compression(),
            mainNtuple->dir().seek_directory(),
            ntupleDescription->GetMainBranches(), basketSizes,
            ntupleDescription->GetDescription().GetNtupleBooking(),
            fRowMode, basketEntries, verbose);

    ntupleDescription->SetNtuple(
      static_cast<tools::wroot::imt_ntuple*>(mtNtuple));
    ntupleDescription->SetBasePNtuple(
      static_cast<tools::wroot::base_pntuple*>(mtNtuple));
  }

  ntupleDescription->GetDescription().SetIsNtupleOwner(true);
  //        // pntuple object is not deleted automatically
  fNtupleVector.push_back(ntupleDescription->GetNtuple());

  Message(kVL3, "create from main", "pntuple", mainNtuple->name());
}

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtupleDescriptionsFromBooking()
{
// Create ntuple descriptions from booking.
// This function is called from the first Fill call.

  // Create pntuple descriptions from ntuple booking.
  auto g4NtupleBookings = fBookingManager->GetNtupleBookingVector();

  for ( auto g4NtupleBooking : g4NtupleBookings ) {
    auto ntupleDescription = new G4RootPNtupleDescription(g4NtupleBooking);
    // Save g4booking, activation in pntuple booking
    fNtupleDescriptionVector.push_back(ntupleDescription);
  }
}

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtuplesFromMain()
{
// Create slave ntuples from  main ntuple.
// This function is called from the first Fill call.

  auto& mainNtupleVector = fMainNtupleManager->GetNtupleVector();

  G4int lcounter = 0;
  for ( auto mainNtuple : mainNtupleVector ) {
    auto& ntupleDescription = fNtupleDescriptionVector[lcounter++];
    CreateNtupleFromMain(ntupleDescription, mainNtuple);
  }
}

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtuplesIfNeeded()
{
// The ntuples on workers are created at first FillColumn or AddRow
// (if only columns of vector type) call.
// When writing multiple times in teh same file, the main ntuples have
// to be recreated as well.

  // G4cout << "G4RootPNtupleManager::CreateNtuplesIfNeeded: "
  //   << " fCreateNtuple: " << fCreateNtuples << " "
  //   << " fNewCycle: " << fNewCycle
  //   << " fMainNtupleManager->GetNewCycle(): " << fMainNtupleManager->GetNewCycle()
  //   << " fNtupleDescriptionVector.size(): " << fNtupleDescriptionVector.size()
  //   << " fNtupleVector.size(): " << fNtupleVector.size()
  //   << G4endl;

  if (fCreateNtuples) {
    // create ntuple descriptions
    CreateNtupleDescriptionsFromBooking();

    // create main ntuples if needed
    G4AutoLock lock(&createMainMutex);
    if (fMainNtupleManager->GetNewCycle()) {
      fMainNtupleManager->CreateNtuplesFromBooking();
    }
    lock.unlock();

    // create slave ntuples
    CreateNtuplesFromMain();
    fCreateNtuples = false;
  }

  if (fNewCycle) {
    // create main ntuples if needed
    G4AutoLock lock(&createMainMutex);
    if (fMainNtupleManager->GetNewCycle()) {
      fMainNtupleManager->CreateNtuplesFromBooking();
    }
    lock.unlock();

    // create slave ntuples
    CreateNtuplesFromMain();
    fNewCycle = false;
  }
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::CreateNtuple(G4NtupleBooking* /*booking*/)
{
// Create pntuple from g4 ntuple booking.
// Nothing to be done here.

  return G4Analysis::kInvalidId;
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

  if ( IsVerbose(kVL4) ) {
    Message(kVL4, "add", "pntuple row", " ntupleId " + to_string(ntupleId));
  }

  // Creating ntuples on workers is triggered with the first FillColumn
  // or AddRow (in only columns of vector type call)
  CreateNtuplesIfNeeded();

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "AddNtupleRow");
  if (ntupleDescription == nullptr) return false;

  auto rfile = std::get<0>(*ntupleDescription->GetDescription().GetFile());

  G4AutoLock lock(&pntupleMutex);
  lock.unlock();
  mutex toolsLock(lock);
  auto result
    = ntupleDescription->GetNtuple()->add_row(toolsLock, *rfile);

  if ( ! result ) {
    Warn("NtupleId " + to_string(ntupleId) + "adding row failed.",
      fkClass, "AddNtupleRow");
  }

  ntupleDescription->GetDescription().SetHasFill(true);

  if ( IsVerbose(kVL3) ) {
    Message(kVL3, "add", "pntuple row", " ntupleId " + to_string(ntupleId));
  }


  return true;
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::Merge()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector) {

    // skip inactivated ntuples
    if (! ntupleDescription->GetDescription().GetActivation() ||
        (ntupleDescription->GetNtuple() == nullptr)) {
      // G4cout << "skipping inactive ntuple " << G4endl;
      continue;
    }

    if ( IsVerbose(kVL4) ) {
      Message(kVL4, "merge", "pntuple",
        ntupleDescription->GetDescription().GetNtupleBooking().name());
    }

    auto rfile = std::get<0>(*ntupleDescription->GetDescription().GetFile());

    G4AutoLock lock(&pntupleMutex);
    lock.unlock();
    mutex toolsLock(lock);
    auto result
      = ntupleDescription->GetNtuple()->end_fill(toolsLock, *rfile);

    if ( ! result ) {
      Warn("Ntuple " + ntupleDescription->GetDescription().GetNtupleBooking().name() +
           "end fill has failed.", fkClass, "Merge");
    }

    if ( IsVerbose(kVL3) ) {
      Message(kVL3, "merge", "pntuple",
        ntupleDescription->GetDescription().GetNtupleBooking().name());
    }

  }

  // Set new cycle
  fNewCycle = true;

  return true;
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::Reset()
{
  // Reset ntuple description, this will delete ntuple if present and
  // we have its ownership.
  // The ntuples will be recreated with new cycle or new open file.

  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    ntupleDescription->Reset();
  }

  fNtupleVector.clear();

  return true;
}

//_____________________________________________________________________________
void G4RootPNtupleManager::Clear()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    delete ntupleDescription->GetNtuple();
  }

  fNtupleDescriptionVector.clear();
  fNtupleVector.clear();

  Message(kVL2, "clear", "pntuples");
}

//_____________________________________________________________________________

void  G4RootPNtupleManager::SetActivation(
  G4bool activation)
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    ntupleDescription->GetDescription().SetActivation(activation);
  }
}

//_____________________________________________________________________________

void  G4RootPNtupleManager::SetActivation(
  G4int ntupleId, G4bool activation)
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetActivation");
  if (ntupleDescription == nullptr) return;

  ntupleDescription->GetDescription().SetActivation(activation);
}

//_____________________________________________________________________________
G4bool  G4RootPNtupleManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "GetActivation");
  if (ntupleDescription == nullptr) return false;

  return ntupleDescription->GetDescription().GetActivation();
}

//_____________________________________________________________________________
void G4RootPNtupleManager::SetNewCycle(G4bool value)
{
  fNewCycle = value;
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::GetNewCycle() const
{
  return fNewCycle;
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::GetNofNtuples() const
{
  return (G4int)fNtupleVector.size();
}

//_____________________________________________________________________________
void G4RootPNtupleManager::SetNtupleRowWise(G4bool rowWise, G4bool rowMode)
{
  fRowWise = rowWise;
  fRowMode = rowMode;
}

//_____________________________________________________________________________
G4bool G4RootPNtupleManager::List(std::ostream& /*output*/, G4bool /*onlyIfActive*/)
{
  Warn("Not implemented.", fkClass, "List");
  return false;
}
