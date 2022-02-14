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

// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#include "G4CsvRNtupleManager.hh"
#include "G4CsvRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4CsvRNtupleManager::G4CsvRNtupleManager(const G4AnalysisManagerState& state)
 : G4VRNtupleManager(state),
   fNtupleVector()
{
}

//
// private methods
//

//_____________________________________________________________________________
G4int G4CsvRNtupleManager::ReadNtupleImpl(const G4String& ntupleName,
                                          const G4String& fileName,
                                          const G4String& dirName,
                                          G4bool isUserFileName)
{
  Message(kVL4, "read", "ntuple", ntupleName);

  // Ntuples are saved per object and per thread
  // but apply the ntuple name and the thread suffixes
  // only if fileName is not provided explicitly
  G4String fullFileName = fileName;
  if ( ! isUserFileName ) {
    fullFileName = fFileManager->GetNtupleFileName(ntupleName);
  }

  // Update directory path
  if ( ! dirName.empty() ) {
    fullFileName = "./" + dirName + "/" + fullFileName;
  }

  // Open file
  if ( ! fFileManager->OpenRFile(fullFileName) ) return kInvalidId;
  auto ntupleFile = fFileManager->GetRFile(fullFileName);

  // Create ntuple
  auto rntuple = new tools::rcsv::ntuple(*ntupleFile);
  auto id = SetNtuple(new G4TRNtupleDescription<tools::rcsv::ntuple>(rntuple));

  Message(kVL2, "read", "ntuple", ntupleName, id > kInvalidId);

  return id;
}

//_____________________________________________________________________________
G4bool G4CsvRNtupleManager::GetTNtupleRow(
  G4TRNtupleDescription<tools::rcsv::ntuple>* ntupleDescription)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL4()->Message("get", "ntuple row", description);
  }  
#endif

  G4CsvRNtupleDescription* ntupleDescription
    = GetNtupleInFunction(ntupleId, "GetNtupleRow");
  if ( ! ntupleDescription )  return false;   
  
  tools::rcsv::ntuple* ntuple = ntupleDescription->fNtuple;

  G4bool isInitialized = ntupleDescription->fIsInitialized;
  if ( ! isInitialized ) {
    tools::ntuple_binding* ntupleBinding = ntupleDescription->fNtupleBinding;
    if ( ! ntuple->initialize(G4cout, *ntupleBinding) ) {
      Warn("Ntuple initialization failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
    ntupleDescription->fIsInitialized = true;
    ntuple->start();
  }

  G4bool next = ntuple->next();
  if ( next ) {
    if ( ! ntuple->get_row() ) {
      Warn("Ntuple get_row() failed !!", fkClass, "GetTNtupleRow");
      return false;
    }
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    G4ExceptionDescription description;
    description << " ntupleId " << ntupleId;  
    fState.GetVerboseL2()->Message("get", "ntuple row", description, true);
  }  
#endif

  return next;
}
