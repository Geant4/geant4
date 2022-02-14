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
#include "G4RootMainNtupleManager.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/ntuple"

using namespace G4Analysis;
using std::to_string;

// mutex in a file scope
namespace {

//Mutex to lock master manager when adding ntuple row and ending fill
G4Mutex pntupleMutex = G4MUTEX_INITIALIZER;

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
  if ( ! ntupleDescription ) return nullptr;

  if ( ! ntupleDescription->fBasePNtuple ) {
    if ( warn ) {
      NotExistWarning("ntuple", id, fkClass, functionName);
    }
    return nullptr;
  }
  return ntupleDescription->fBasePNtuple;
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

  auto file = fMainNtupleManager->GetNtupleFile(&ntupleDescription->fDescription);
  if ( ! file ) {
    Warn("Cannot create pntuple. Main ntuple file does not exist.",
      fkClass, "CreateNtupleFromMain");
    return;
  }

  ntupleDescription->fDescription.fFile = file;

  // Get parameters from ntupleDescription
  mainNtuple->get_branches(ntupleDescription->fMainBranches);

  auto rfile = std::get<0>(*file);
  G4bool verbose = true;
  if ( fRowWise ) {
    auto mainBranch = mainNtuple->get_row_wise_branch();
    auto mtNtuple
      = new tools::wroot::mt_ntuple_row_wise(
              G4cout, rfile->byte_swap(), rfile->compression(),
              mainNtuple->dir().seek_directory(),
              *mainBranch, mainBranch->basket_size(),
              ntupleDescription->fDescription.fNtupleBooking, verbose);

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
    auto basketEntries = fMainNtupleManager->GetBasketEntries();

    auto mtNtuple =
      new tools::wroot::mt_ntuple_column_wise(
            G4cout, rfile->byte_swap(), rfile->compression(),
            mainNtuple->dir().seek_directory(),
            ntupleDescription->fMainBranches, basketSizes,
            ntupleDescription->fDescription.fNtupleBooking,
            fRowMode, basketEntries, verbose);

    ntupleDescription->fNtuple
      = static_cast<tools::wroot::imt_ntuple*>(mtNtuple);
    ntupleDescription->fBasePNtuple
      = static_cast<tools::wroot::base_pntuple*>(mtNtuple);
  }

  ntupleDescription->fDescription.fIsNtupleOwner = true;
  //        // pntuple object is not deleted automatically
  fNtupleVector.push_back(ntupleDescription->fNtuple);

  Message(kVL3, "create from main", "pntuple", mainNtuple->name());
}

//_____________________________________________________________________________
void G4RootPNtupleManager::CreateNtuplesFromMain()
{
// Create ntuple from booking (if not yet done) and main ntuple
// This function is called from the first Fill call.

  // Create pntuple descriptions from ntuple booking.
  auto g4NtupleBookings = fBookingManager->GetNtupleBookingVector();
  for ( auto g4NtupleBooking : g4NtupleBookings ) {
    auto ntupleDescription = new G4RootPNtupleDescription(g4NtupleBooking);
    // Save g4booking, activation in pntuple booking
    fNtupleDescriptionVector.push_back(ntupleDescription);
  }

  auto& mainNtupleVector = fMainNtupleManager->GetNtupleVector();

  G4int lcounter = 0;
  for ( auto mainNtuple : mainNtupleVector ) {
    auto& ntupleDescription = fNtupleDescriptionVector[lcounter++];
    CreateNtupleFromMain(ntupleDescription, mainNtuple);
  }

  fCreateNtuples = false;
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
  if (fCreateNtuples) {
    CreateNtuplesFromMain();
  }

  if ( fState.GetIsActivation() && ( ! GetActivation(ntupleId) ) ) {
    //G4cout << "Skipping AddNtupleRow for " << ntupleId << G4endl;
    return false;
  }

  if ( IsVerbose(kVL4) ) {
    Message(kVL4, "add", "pntuple row", " ntupleId " + to_string(ntupleId));
  }

  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "AddNtupleRow");
  if ( ! ntupleDescription ) return false;

  auto rfile = std::get<0>(*ntupleDescription->fDescription.fFile);

  G4AutoLock lock(&pntupleMutex);
  lock.unlock();
  mutex toolsLock(lock);
  auto result
    = ntupleDescription->fNtuple->add_row(toolsLock, *rfile);

  if ( ! result ) {
    Warn("NtupleId " + to_string(ntupleId) + "adding row failed.",
      fkClass, "AddNtupleRow");
  }

  ntupleDescription->fDescription.fHasFill = true;

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
    if(!ntupleDescription->fDescription.fActivation || !ntupleDescription->fNtuple) {
      // G4cout << "skipping inactive ntuple " << G4endl;
      continue;
    }

    if ( IsVerbose(kVL4) ) {
      Message(kVL4, "merge", "pntuple", ntupleDescription->fDescription.fNtupleBooking.name());
    }

    auto rfile = std::get<0>(*ntupleDescription->fDescription.fFile);

    G4AutoLock lock(&pntupleMutex);
    lock.unlock();
    mutex toolsLock(lock);
    auto result
      = ntupleDescription->fNtuple->end_fill(toolsLock, *rfile);

    if ( ! result ) {
      Warn("Ntuple " + ntupleDescription->fDescription.fNtupleBooking.name() +
           "end fill has failed.", fkClass, "Merge");
    }

    delete ntupleDescription->fNtuple;
    ntupleDescription->fNtuple = nullptr;

    if ( IsVerbose(kVL3) ) {
      Message(kVL3, "merge", "pntuple", ntupleDescription->fDescription.fNtupleBooking.name());
    }

  }
  return true;

}

//_____________________________________________________________________________
void G4RootPNtupleManager::Clear()
{
  for ( auto ntupleDescription : fNtupleDescriptionVector ) {
    delete ntupleDescription->fNtuple;
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
    ntupleDescription->fDescription.fActivation = activation;
  }
}

//_____________________________________________________________________________

void  G4RootPNtupleManager::SetActivation(
  G4int ntupleId, G4bool activation)
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "SetActivation");
  if ( ! ntupleDescription ) return;

  ntupleDescription->fDescription.fActivation = activation;
}

//_____________________________________________________________________________
G4bool  G4RootPNtupleManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleDescription = GetNtupleDescriptionInFunction(ntupleId, "GetActivation");
  if ( ! ntupleDescription ) return false;

  return ntupleDescription->fDescription.fActivation;
}

//_____________________________________________________________________________
G4int G4RootPNtupleManager::GetNofNtuples() const
{
  return fNtupleVector.size();
}

//_____________________________________________________________________________
void G4RootPNtupleManager::SetNtupleRowWise(G4bool rowWise, G4bool rowMode)
{
  fRowWise = rowWise;
  fRowMode = rowMode;
}
