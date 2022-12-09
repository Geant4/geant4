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

// Author: Ivana Hrivnacova, 01/09/2020  (ivana@ipno.in2p3.fr)

#include "G4NtupleBookingManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;
using std::to_string;

//
// private template functions
//

//_____________________________________________________________________________
G4NtupleBookingManager::G4NtupleBookingManager(
  const G4AnalysisManagerState& state)
  : G4BaseAnalysisManager(state)
{}

//_____________________________________________________________________________
G4NtupleBookingManager::~G4NtupleBookingManager()
{
  for ( auto ntupleBooking : fNtupleBookingVector ) {
    delete ntupleBooking;
  }
}

//_____________________________________________________________________________
G4NtupleBooking*
G4NtupleBookingManager::GetNtupleBookingInFunction(
  G4int id, std::string_view functionName, G4bool warn) const
{
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleBookingVector.size()) ) {
    if ( warn) {
      Warn("Ntuple booking " + to_string(id) + " does not exist.",
        fkClass, functionName);
    }
    return nullptr;
  }

  return fNtupleBookingVector[index];
}

//_____________________________________________________________________________
G4bool G4NtupleBookingManager::CheckName(
  const G4String& name, const G4String& objectType) const
{
  if (name.size() == 0u) {
    Warn("Empty " + objectType + " name is not allowed.\n" +
      objectType + " was not created.",  fkClass, "CheckName");
    return false;
  }
  return true;
}

//
// protected functions
//

//_____________________________________________________________________________
G4bool G4NtupleBookingManager::IsEmpty() const
{
  return fNtupleBookingVector.size() == 0u;
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtuple(
  const G4String& name, const G4String& title)
{
  if ( ! CheckName(name, "Ntuple") ) return kInvalidId;

  Message(kVL4, "create", "ntuple booking", name);

  // Create ntuple description
  auto index = fNtupleBookingVector.size();
  auto ntupleBooking = new G4NtupleBooking();
  fNtupleBookingVector.push_back(ntupleBooking);

  // Save name & title in ntuple booking
  ntupleBooking->fNtupleBooking.set_name(name);
  ntupleBooking->fNtupleBooking.set_title(title);
  ntupleBooking->fNtupleId = G4int(index + fFirstId);

  // Create ntuple
  fLockFirstId = true;

  Message(kVL2, "create", "ntuple booking",
    name +  " ntupleId " + to_string(ntupleBooking->fNtupleId));

  return ntupleBooking->fNtupleId;
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleIColumn(const G4String& name,
                                               std::vector<int>* vector)
{
  return CreateNtupleIColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleFColumn(const G4String& name,
                                               std::vector<float>* vector)
{
  return CreateNtupleFColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleDColumn(const G4String& name,
                                               std::vector<double>* vector)
{
  return CreateNtupleDColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleSColumn(const G4String& name,
                                               std::vector<std::string>* vector)
{
  return CreateNtupleSColumn(GetCurrentNtupleId(), name, vector);
}

//_____________________________________________________________________________
G4NtupleBooking*  G4NtupleBookingManager::FinishNtuple()
{
  return FinishNtuple(GetCurrentNtupleId());
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleIColumn(
  G4int ntupleId, const G4String& name, std::vector<int>* vector)
{
  return CreateNtupleTColumn<int>(ntupleId, name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleFColumn(
  G4int ntupleId, const G4String& name, std::vector<float>* vector)
{
  return CreateNtupleTColumn<float>(ntupleId, name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleDColumn(
  G4int ntupleId, const G4String& name, std::vector<double>* vector)
{
  return CreateNtupleTColumn<double>(ntupleId, name, vector);
}

//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtupleSColumn(
  G4int ntupleId, const G4String& name, std::vector<std::string>* vector)
{
  return CreateNtupleTColumn<std::string>(ntupleId, name, vector);
}

//_____________________________________________________________________________
G4NtupleBooking*  G4NtupleBookingManager::FinishNtuple(
  G4int ntupleId)
{
  // Nothing to be done for booking,
  return GetNtupleBookingInFunction(ntupleId, "FinishNtuple");
}

//_____________________________________________________________________________
G4bool G4NtupleBookingManager::SetFirstNtupleColumnId(G4int firstId)
{
  if ( fLockFirstNtupleColumnId ) {
    Warn("Cannot set FirstNtupleColumnId as its value was already used.",
      fkClass, "SetFirstNtupleColumnId");
    return false;
  }

  fFirstNtupleColumnId = firstId;
  return true;
}

//_____________________________________________________________________________
void  G4NtupleBookingManager::SetActivation(
  G4bool activation)
{
  for ( auto ntupleBooking : fNtupleBookingVector ) {
    ntupleBooking->fActivation = activation;
  }
}

//_____________________________________________________________________________
void  G4NtupleBookingManager::SetActivation(
  G4int ntupleId, G4bool activation)
{
  auto ntupleBooking
    = GetNtupleBookingInFunction(ntupleId, "SetActivation");
  if (ntupleBooking == nullptr) return;

  ntupleBooking->fActivation = activation;
}

//_____________________________________________________________________________
G4bool  G4NtupleBookingManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleBooking
    = GetNtupleBookingInFunction(ntupleId, "GetActivation");
  if (ntupleBooking == nullptr) return false;

  return ntupleBooking->fActivation;
}

//_____________________________________________________________________________
void  G4NtupleBookingManager::SetFileName(
  const G4String& fileName)
{
  for ( auto ntupleBooking : fNtupleBookingVector ) {
    ntupleBooking->fFileName = fileName;
  }
}

//_____________________________________________________________________________
void  G4NtupleBookingManager::SetFileName(
  G4int ntupleId, const G4String& fileName)
{
  auto ntupleBooking
    = GetNtupleBookingInFunction(ntupleId, "SetFileName");
  if (ntupleBooking == nullptr) return;

  // Do nothing if file name does not change
  if ( ntupleBooking->fFileName == fileName ) return;

  auto ntupleFileName = fileName;
  auto extension = GetExtension(fileName);
  if (extension.size() != 0u) {
    // Check if valid extension (if present)
    auto output = G4Analysis::GetOutput(extension);
    if ( output == G4AnalysisOutput::kNone ) {
      Warn("The file extension " + extension + " is not supported.",
        fkClass, "SetFileName");
      return;
    }
  }
  else {
    if (fFileType.size() != 0u) {
      //add extension if missing and file type is defined
      ntupleFileName = fileName + "." + fFileType;
    }
  }

  // Save the fileName in booking
  // If extension is still missing (possible with generic manager),
  // it will be completed with the default one at OpenFile
  ntupleBooking->fFileName = ntupleFileName;
}

//_____________________________________________________________________________
G4String  G4NtupleBookingManager::GetFileName(
  G4int ntupleId) const
{
  auto ntupleBooking
    = GetNtupleBookingInFunction(ntupleId, "GetFileName");
  if (ntupleBooking == nullptr) return "";

  return ntupleBooking->fFileName;
}

//_____________________________________________________________________________
void G4NtupleBookingManager::ClearData()
{
  for ( auto ntupleBooking : fNtupleBookingVector ) {
    delete ntupleBooking;
  }
  fNtupleBookingVector.clear();
  fLockFirstNtupleColumnId = false;

  Message(G4Analysis::kVL2, "clear", "ntupleBookings");
}

//
// public methods
//

//_____________________________________________________________________________
void G4NtupleBookingManager::SetFileType(const G4String& fileType)
{
  // do nothing if file type is defined and is same
  if ( fFileType == fileType ) return;

  // save the type
  fFileType = fileType;

  // Give warning and redefine file extension in bookings
  // with file name of different fileTypes
  for ( auto ntupleBooking : fNtupleBookingVector ) {
    if ((ntupleBooking->fFileName).size() == 0u) continue;

    auto extension = GetExtension(ntupleBooking->fFileName);
    if ( fFileType == extension ) continue;

    // multiple file types are not suported
    auto baseFileName = GetBaseName(ntupleBooking->fFileName);
    auto ntupleFileName = baseFileName + "." + fFileType;
    if (extension.size() != 0u) {
      Warn("Writing ntuples in files of different output types " +
           fFileType + ", " + extension + " is not supported.",
           fkClass, "SetFileType");
    }

    // Save the info in ntuple description
    ntupleBooking->fFileName = ntupleFileName;
  }
}
