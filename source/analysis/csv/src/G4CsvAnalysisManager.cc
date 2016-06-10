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
// $Id: G4CsvAnalysisManager.cc 74257 2013-10-02 14:24:55Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisManager.hh"
#include "G4CsvFileManager.hh"
#include "G4H1DummyManager.hh"
#include "G4H2DummyManager.hh"
#include "G4CsvNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>

G4CsvAnalysisManager* G4CsvAnalysisManager::fgMasterInstance = 0;
G4ThreadLocal G4CsvAnalysisManager* G4CsvAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4CsvAnalysisManager* G4CsvAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4CsvAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager(G4bool isMaster)
 : G4VAnalysisManager("Csv", isMaster),
   fNtupleManager(0),
   fFileManager(0)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description << "      " 
                << "G4CsvAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("G4CsvAnalysisManager::G4CsvAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;
  
  // Create managers
  fH1Manager = new G4H1DummyManager(fState);
  fH2Manager = new G4H2DummyManager(fState);
  fNtupleManager = new G4CsvNtupleManager(fState);
  fFileManager = new G4CsvFileManager(fState);
  fNtupleManager->SetFileManager(fFileManager);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{  
  if ( fState.GetIsMaster() ) fgMasterInstance = 0;
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseNtupleFiles()
{
  const std::vector<G4CsvNtupleDescription*>& ntupleVector
    = fNtupleManager->GetNtupleVector();

  // Close ntuple files
  std::vector<G4CsvNtupleDescription*>::const_iterator it;  
  for (it = ntupleVector.begin(); it != ntupleVector.end(); it++ ) {
    fFileManager->CloseNtupleFile((*it));
  }
  
  return true;
}    


// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;

  // Only book file name in file manager
  G4bool result = fFileManager->OpenFile(fileName);
  finalResult = finalResult && result;

  // Create ntuples if they are booked  
  // (The files will be created with creating ntuples)
  fNtupleManager->CreateNtuplesFromBooking();


  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteImpl() 
{
  // nothing to be done for Csv file
  G4bool result = true;
  
  // Write ASCII if activated
  // Not available 
  //if ( IsAscii() ) {
  //  result = WriteAscii();
  //}   

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFileImpl()
{
  G4bool finalResult = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "files", "");
#endif

  // Unlock file name only
  G4bool result = fFileManager->CloseFile();
  finalResult = finalResult && result;
  
  // Close ntuple files
  if ( ( ! G4AnalysisManagerState::IsMT() ) || ( ! fState.GetIsMaster() ) ) {
    // In sequential mode or in MT mode only on workers
    result = CloseNtupleFiles();
    finalResult = finalResult && result;
  }  

  // reset data
  result = fNtupleManager->Reset();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Resetting data failed";
    G4Exception("G4CsvAnalysisManager::CloseFile()",
              "Analysis_W002", JustWarning, description);
    result = false;       
  } 
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) {
    fState.GetVerboseL1()->Message("close", "files", "");
  }  
#endif

  return finalResult; 
} 
   
