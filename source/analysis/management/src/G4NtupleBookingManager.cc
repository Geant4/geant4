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

//
// private template functions
//

//_____________________________________________________________________________
G4NtupleBookingManager::G4NtupleBookingManager(
  const G4AnalysisManagerState& state)
  : G4BaseAnalysisManager(state),
    fNtupleBookingVector(),
    fFileType(),
    fFirstNtupleColumnId(0),
    fLockFirstNtupleColumnId(false)    
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
  G4int id, G4String functionName, G4bool warn) const
{                                      
  auto index = id - fFirstId;
  if ( index < 0 || index >= G4int(fNtupleBookingVector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4NtupleBookingManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "ntuple booking " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return nullptr;         
  }
  
  return fNtupleBookingVector[index];
}

//
// protected functions
//

//_____________________________________________________________________________
G4bool G4NtupleBookingManager::IsEmpty() const
{
  return ! fNtupleBookingVector.size();
}  
 
//_____________________________________________________________________________
G4int G4NtupleBookingManager::CreateNtuple(
  const G4String& name, const G4String& title)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "ntuple booking", name);
#endif

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

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << name << " ntupleId " << ntupleBooking->fNtupleId;
    fState.GetVerboseL2()->Message("create", "ntuple booking", description);
  } 
#endif

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
G4int G4NtupleBookingManager::CreateNtupleSColumn(const G4String& name)
{
  return CreateNtupleSColumn(GetCurrentNtupleId(), name);
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
  G4int ntupleId, const G4String& name)
{
  return CreateNtupleTColumn<std::string>(ntupleId, name, nullptr);
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
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstNtupleColumnId as its value was already used.";
    G4Exception("G4BaseNtupleManager::SetFirstNtupleColumnId()",
                "Analysis_W013", JustWarning, description);
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
  if ( ! ntupleBooking ) return;

  ntupleBooking->fActivation = activation;
}

//_____________________________________________________________________________
G4bool  G4NtupleBookingManager::GetActivation(
  G4int ntupleId) const
{
  auto ntupleBooking 
    = GetNtupleBookingInFunction(ntupleId, "GetActivation");
  if ( ! ntupleBooking ) return false;

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
  if ( ! ntupleBooking ) return;

  // Do nothing if file name does not change
  if ( ntupleBooking->fFileName == fileName ) return;

  auto ntupleFileName = fileName;
  auto extension = GetExtension(fileName);
  if ( extension.size() ) {
    // Check if valid extension (if present)
    auto output = G4Analysis::GetOutput(extension);
    if ( output == G4AnalysisOutput::kNone ) {
      G4ExceptionDescription description;
      description
        << "The file extension " << extension << "is not supported.";
      G4Exception("G4NtupleBookingManager::SetFileName",
                  "Analysis_W051", JustWarning, description);
      return;
    }
  }
  else {
    if ( fFileType.size() ) {
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
  if ( ! ntupleBooking ) return "";

  return ntupleBooking->fFileName;
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

    if ( ! (ntupleBooking->fFileName).size() ) continue;

    auto extension = GetExtension(ntupleBooking->fFileName);
    if ( fFileType == extension ) continue;

    // multiple file types are not suported
    auto baseFileName = GetBaseName(ntupleBooking->fFileName);
    auto ntupleFileName = baseFileName + "." + fFileType;
    if ( extension.size()) {
    G4ExceptionDescription description;
      description
        << "Writing ntuples in files of different output types "
        << fFileType << ", " << extension << " is not supported." << G4endl
        << "Ntuple " << ntupleBooking->fNtupleBooking.name()
        << " will be written in " << ntupleFileName;
      G4Exception("G4NtupleBookingManager::SetFileType",
                  "Analysis_W051", JustWarning, description);
    }

    // Save the info in ntuple description
    ntupleBooking->fFileName = ntupleFileName;
  }
}
