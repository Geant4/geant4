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

#include "G4CsvAnalysisManager.hh"
#include "G4CsvFileManager.hh"
#include "G4CsvNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH2Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH3Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging P1 profiles
  G4Mutex mergeP1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging P2 profiles
  G4Mutex mergeP2Mutex = G4MUTEX_INITIALIZER;
}  

G4CsvAnalysisManager* G4CsvAnalysisManager::fgMasterInstance = nullptr;
G4ThreadLocal G4CsvAnalysisManager* G4CsvAnalysisManager::fgInstance = nullptr;

//_____________________________________________________________________________
G4CsvAnalysisManager* G4CsvAnalysisManager::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4CsvAnalysisManager(isMaster);
  }
  
  return fgInstance;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::IsInstance()
{
  return ( fgInstance != 0 );
}    

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager(G4bool isMaster)
 : G4ToolsAnalysisManager("Csv", isMaster),
   fNtupleManager(nullptr),
   fFileManager(nullptr)
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
  fNtupleManager = new G4CsvNtupleManager(fState);
  fFileManager = std::make_shared<G4CsvFileManager>(fState);
  fNtupleManager->SetFileManager(fFileManager);
      // The managers will be deleted by the base class
  
  // Set managers to base class which takes then their ownreship
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{  
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseNtupleFiles()
{
  auto ntupleVector
    = fNtupleManager->GetNtupleDescriptionVector();

  // Close ntuple files
  for ( auto ntupleDescription : ntupleVector) {
    fFileManager->CloseNtupleFile(ntupleDescription);
  }

  return true;
}    


// 
// protected methods
//

// 
// private methods
//


//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteH1()
{
  auto h1Vector = fH1Manager->GetH1Vector();
  auto hnVector = fH1Manager->GetHnVector();

  if ( ! h1Vector.size() ) return true;

  auto result = true;

  if ( ! G4Threading::IsWorkerThread() )  {
    result = WriteT(h1Vector, hnVector, "h1");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH1(&mergeH1Mutex);
    fgMasterInstance->fH1Manager->AddH1Vector(h1Vector);
    lH1.unlock();
  }  
  
  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteH2()
{
  auto h2Vector = fH2Manager->GetH2Vector();
  auto hnVector = fH2Manager->GetHnVector();

  if ( ! h2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    result = WriteT(h2Vector, hnVector, "h2");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH2(&mergeH2Mutex);
    fgMasterInstance->fH2Manager->AddH2Vector(h2Vector);
    lH2.unlock();
  }  
  
  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteH3()
{
  auto h3Vector = fH3Manager->GetH3Vector();
  auto hnVector = fH3Manager->GetHnVector();

  if ( ! h3Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    result = WriteT(h3Vector, hnVector, "h3");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH3(&mergeH3Mutex);
    fgMasterInstance->fH3Manager->AddH3Vector(h3Vector);
    lH3.unlock();
  }  
  
  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteP1()
{
  auto p1Vector = fP1Manager->GetP1Vector();
  auto hnVector = fP1Manager->GetHnVector();

  if ( ! p1Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    result = WriteT(p1Vector, hnVector, "p1");
  }  
  else {
    // The worker manager just adds its profiles to the master
    // This operation needs a lock
    G4AutoLock lP1(&mergeP1Mutex);
    fgMasterInstance->fP1Manager->AddP1Vector(p1Vector);
    lP1.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteP2()
{
  auto p2Vector = fP2Manager->GetP2Vector();
  auto hnVector = fP2Manager->GetHnVector();

  if ( ! p2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    result = WriteT(p2Vector, hnVector, "p2");
  }  
  else {
    // The worker manager just adds its profiles to the master
    // This operation needs a lock
    G4AutoLock lP2(&mergeP2Mutex);
    fgMasterInstance->fP2Manager->AddP2Vector(p2Vector);
    lP2.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;

  auto result = G4ToolsAnalysisManager::Reset();
  finalResult = finalResult && result;

  result = fNtupleManager->Reset(true);
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  auto finalResult = true;
  auto result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

  // Only lock file name in file manager
  result = fFileManager->OpenFile(fileName);
  finalResult = finalResult && result;

  // Histogram and profile files are created/closed indivudually for each
  // object in WriteHn{Pn] function

  // Create ntuples if they are booked  
  // (The files will be created with creating ntuples)
  fNtupleManager->CreateNtuplesFromBooking();

  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteImpl() 
{
  // nothing to be done for Csv file
  auto finalResult = true;
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write", "files", "");
#endif


  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) || 
         ( ! fH3Manager->IsEmpty() ) || ( ! fP1Manager->IsEmpty() ) ||
         ( ! fP2Manager->IsEmpty() ) ) ) {

    G4ExceptionDescription description;
    description 
      << "      " << "No master G4CsvAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram data will not be merged.";
      G4Exception("G4CsvAnalysisManager::Write()",
                "Analysis_W031", JustWarning, description);
  }

  // H1
  auto result = WriteH1();
  finalResult = finalResult && result;

  // H2
  result = WriteH2();
  finalResult = finalResult && result;

  // H3
  result = WriteH3();
  finalResult = finalResult && result;

  // P1
  result = WriteP1();
  finalResult = finalResult && result;

  // P2
  result = WriteP2();
  finalResult = finalResult && result;
  
  // Ntuples
  // Nothing to be done

  // Write ASCII if activated
  // Not available 
  //if ( IsAscii() ) {
  //  result = WriteAscii();
  //}   

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("write", "files", "", finalResult);
#endif

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFileImpl(G4bool reset)
{
  auto finalResult = true;

  // Unlock file name only
  auto result = fFileManager->CloseFile();
  finalResult = finalResult && result;
   
  // Histogram and profile files are created/closed indivudually for each
  // object in WriteHn{Pn] function
  // In sequential mode or in MT mode only on workers
  result = CloseNtupleFiles();
  finalResult = finalResult && result;

  // reset data
  if ( reset ) {
    result = Reset();
  } else {
    // ntuple must be reset 
    result = fNtupleManager->Reset(true);
  }
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Resetting data failed";
    G4Exception("G4CsvAnalysisManager::CloseFile()",
              "Analysis_W021", JustWarning, description);
  } 
  finalResult = finalResult && result;

  return finalResult; 
} 
