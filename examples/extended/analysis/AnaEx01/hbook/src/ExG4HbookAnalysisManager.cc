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
/// \file hbook/src/ExG4HbookAnalysisManager.cc
/// \brief Implementation of the ExG4HbookAnalysisManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookAnalysisManager.hh"
#include "ExG4HbookFileManager.hh"
#include "ExG4HbookH1Manager.hh"
#include "ExG4HbookH2Manager.hh"
#include "ExG4HbookH3DummyManager.hh"
#include "ExG4HbookP1Manager.hh"
#include "ExG4HbookP2DummyManager.hh"
#include "ExG4HbookNtupleManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4UnitsTable.hh"

#include <iostream>

extern "C" int setpawc();
extern "C" int setntuc();

ExG4HbookAnalysisManager* ExG4HbookAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
ExG4HbookAnalysisManager* ExG4HbookAnalysisManager::Create(G4bool /*isMaster*/)
{
  if ( fgInstance == 0 ) {
    fgInstance = new ExG4HbookAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
ExG4HbookAnalysisManager* ExG4HbookAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new ExG4HbookAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
ExG4HbookAnalysisManager::ExG4HbookAnalysisManager()
 : G4VAnalysisManager("Hbook", true),
   fH1Manager(0),
   fH2Manager(0),
   fH3Manager(0),
   fP1Manager(0),
   fP2Manager(0),
   fNtupleManager(0),
   fFileManager(0)
{
  if ( G4Threading::IsWorkerThread() ) {
    // Hbook output is not supported in MT mode
    G4ExceptionDescription description;
    description << "      " 
      << "G4HbookAnalysisManager is not supported in multi-threading mode."; 
    G4Exception("ExG4HbookAnalysisManager::ExG4HbookAnalysisManager()",
                "Analysis_F002", FatalException, description);
  }              

  if ( fgInstance ) {
    G4ExceptionDescription description;
    description << "      " 
                << "G4HbookAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("ExG4HbookAnalysisManager::ExG4HbookAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  fgInstance = this;

  // Create managers
  fH1Manager = new ExG4HbookH1Manager(fState);
  fH2Manager = new ExG4HbookH2Manager(fState);
  fH3Manager = new ExG4HbookH3DummyManager(fState);
  fP1Manager = new ExG4HbookP1Manager(fState);
  fP2Manager = new ExG4HbookP2DummyManager(fState);
  fNtupleManager = new ExG4HbookNtupleManager(fState);
  fFileManager = new ExG4HbookFileManager(fState);
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetH3Manager(fH3Manager);
  SetP1Manager(fP1Manager);
  SetP2Manager(fP2Manager);
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
  
  // Set file manager to component managers
  fH1Manager->SetFileManager(fFileManager);
  fH2Manager->SetFileManager(fFileManager);
  fP1Manager->SetFileManager(fFileManager);
  fNtupleManager->SetFileManager(fFileManager);
  
  // Initialize HBOOK :
  tools::hbook::CHLIMIT(setpawc());
  setntuc(); //for ntuple.
}

//_____________________________________________________________________________
ExG4HbookAnalysisManager::~ExG4HbookAnalysisManager()
{  
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::Reset()
{
// Reset histograms and ntuple  

  fH1Manager->Reset();
  fH2Manager->Reset();
  fP1Manager->Reset();
  fNtupleManager->Reset();
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;
  G4bool result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  G4String name = fFileManager->GetFullFileName();
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "analysis file", name);
#endif

  // Open file
  result = fFileManager->OpenFile(fileName);
  finalResult = finalResult && result;

  // Create h1 histrograms 
  fH1Manager->CreateH1sFromBooking();

  // Create h2 histrograms if any is booked
  fH2Manager->CreateH2sFromBooking();

  // Create p1 profiles if any is booked
  fP1Manager->CreateP1sFromBooking();

  // Create ntuples if they are booked
  fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("open", "analysis file", name);
#endif
  
  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::WriteImpl() 
{
  G4bool finalResult = true;

  G4bool result = fFileManager->WriteFile();
  finalResult = finalResult && result;

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii(fFileManager->GetFileName());
    finalResult = finalResult && result;
  }   

  return finalResult;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::CloseFileImpl()
{
  G4bool finalResult = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "file", fFileManager->GetFullFileName());
#endif

  // reset data
  Reset();

  // close file
  G4bool result = fFileManager->CloseFile();
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("close", "file", fFileManager->GetFullFileName(), result);
#endif

  return finalResult;
} 
   
#endif
