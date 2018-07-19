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
// $Id: G4CsvFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvFileManager.hh"
#include "G4AnalysisManagerState.hh"

//_____________________________________________________________________________
G4CsvFileManager::G4CsvFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state)
{}

//_____________________________________________________________________________
G4CsvFileManager::~G4CsvFileManager()
{}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvFileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;

  fLockFileName = true;
  fLockNtupleDirectoryName = true;
  fIsOpenFile = true;
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4CsvFileManager::WriteFile() 
{
  // nothing to be done for Csv file
  return true;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::CloseFile()
{
  fLockFileName = false;
  fIsOpenFile = false;
  return true; 
} 
   
//_____________________________________________________________________________
G4bool G4CsvFileManager::CreateNtupleFile(
  G4TNtupleDescription<tools::wcsv::ntuple>* ntupleDescription)
{
  G4String ntupleName = ntupleDescription->fNtupleBooking.name();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "file", GetNtupleFileName(ntupleName));
#endif

  auto ntupleFile = new std::ofstream(GetNtupleFileName(ntupleName));
  if ( ntupleFile->fail() ) {
    delete ntupleFile;
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " 
                << GetNtupleFileName(ntupleName);
    G4Exception("G4CsvFileManager::CreateNtupleFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("create", "file", GetNtupleFileName(ntupleName));
#endif

  ntupleDescription->fFile = ntupleFile;
  return true;
}  

//_____________________________________________________________________________
G4bool G4CsvFileManager::CloseNtupleFile(
  G4TNtupleDescription<tools::wcsv::ntuple>* ntupleDescription)
{
  // Do nothing if there is no file
  if ( ! ntupleDescription->fFile ) return true;

  G4String ntupleName = ntupleDescription->fNtupleBooking.name();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("close", "file", GetNtupleFileName(ntupleName));
#endif

  // close file
  ntupleDescription->fFile->close(); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("close", "file", GetNtupleFileName(ntupleName));
#endif

  return true; 
} 
   
