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
// $Id: G4RootAnalysisManager.cc 93815 2015-11-02 11:32:20Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootAnalysisManager.hh"
#include "G4RootFileManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <iostream>
#include <cstdio>

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

G4RootAnalysisManager* G4RootAnalysisManager::fgMasterInstance = nullptr;
G4ThreadLocal G4RootAnalysisManager* G4RootAnalysisManager::fgInstance = nullptr;

//_____________________________________________________________________________
G4RootAnalysisManager* G4RootAnalysisManager::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4RootAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::IsInstance()
{
  return ( fgInstance != 0 );
}    

//_____________________________________________________________________________
G4RootAnalysisManager::G4RootAnalysisManager(G4bool isMaster)
 : G4ToolsAnalysisManager("Root", isMaster),
   fNtupleManager(nullptr),
   fFileManager(nullptr)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4RootAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4RootAnalysisManager::G4RootAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // Create managers
  fNtupleManager = new G4RootNtupleManager(fState);
  fFileManager = std::make_shared<G4RootFileManager>(fState);
      // The managers will be deleted by the base class
  
  // Set managers to base class which takes then their ownership
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4RootAnalysisManager::~G4RootAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH1()
{
  auto h1Vector = fH1Manager->GetH1Vector();
  auto hnVector = fH1Manager->GetHnVector();

  if ( ! h1Vector.size() ) return true;

  auto result = true;

  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h1Vector, hnVector, directory, "h1");
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
G4bool G4RootAnalysisManager::WriteH2()
{
  auto h2Vector = fH2Manager->GetH2Vector();
  auto hnVector = fH2Manager->GetHnVector();

  if ( ! h2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h2Vector, hnVector, directory, "h2");
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
G4bool G4RootAnalysisManager::WriteH3()
{
  auto h3Vector = fH3Manager->GetH3Vector();
  auto hnVector = fH3Manager->GetHnVector();

  if ( ! h3Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(h3Vector, hnVector, directory, "h3");
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
G4bool G4RootAnalysisManager::WriteP1()
{
  auto p1Vector = fP1Manager->GetP1Vector();
  auto hnVector = fP1Manager->GetHnVector();

  if ( ! p1Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(p1Vector, hnVector, directory, "p1");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lP1(&mergeP1Mutex);
    fgMasterInstance->fP1Manager->AddP1Vector(p1Vector);
    lP1.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteP2()
{
  auto p2Vector = fP2Manager->GetP2Vector();
  auto hnVector = fP2Manager->GetHnVector();

  if ( ! p2Vector.size() ) return true;

  auto result = true;
  
  if ( ! G4Threading::IsWorkerThread() )  {
    auto directory = fFileManager->GetHistoDirectory(); 
    result = WriteT(p2Vector, hnVector, directory, "p2");
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lP2(&mergeP2Mutex);
    fgMasterInstance->fP2Manager->AddP2Vector(p2Vector);
    lP2.unlock();
  }  
  
  return result;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;
  
  auto result = G4ToolsAnalysisManager::Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset(false);
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  auto finalResult = true;
  auto result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  G4String name = fFileManager->GetFullFileName();
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "analysis file", name);
#endif

  result = fFileManager->OpenFile(fileName);
  finalResult = finalResult && result;
  
  fNtupleManager->SetNtupleDirectory(fFileManager->GetNtupleDirectory());
  fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("open", "analysis file", name);
#endif
  
  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteImpl() 
{

  auto finalResult = true;

  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) || 
         ( ! fH3Manager->IsEmpty() ) || ( ! fP1Manager->IsEmpty() ) || 
         ( ! fP2Manager->IsEmpty() ) ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "No master G4RootAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram/profile data will not be merged.";
      G4Exception("G4RootAnalysisManager::Write()",
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

  // File
  result = fFileManager->WriteFile();
  finalResult = finalResult && result;

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii(fFileManager->GetFileName());
    finalResult = finalResult && result;
  }   

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::CloseFileImpl()
{
  auto finalResult = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("close", "file", fFileManager->GetFullFileName());
#endif

  // reset data
  auto result = Reset();
  if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W021", JustWarning, description);
  } 
  finalResult = finalResult && result;

  // close file
  fFileManager->CloseFile();  

  // No files clean-up in sequential mode
  if ( ! G4Threading::IsMultithreadedApplication() )  return finalResult;
  
  // Delete files if empty in MT mode
  if ( ( fState.GetIsMaster() && 
         fH1Manager->IsEmpty() && fH2Manager->IsEmpty() && fH3Manager->IsEmpty() &&
         fP1Manager->IsEmpty() && fP2Manager->IsEmpty() && fNtupleManager->IsEmpty() ) || 
       ( ( ! fState.GetIsMaster() ) && fNtupleManager->IsEmpty() ) ) {
    result = ! std::remove(fFileManager->GetFullFileName());
    //  std::remove returns 0 when success
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Removing file " 
                  << fFileManager->GetFullFileName() << " failed";
      G4Exception("G4XmlAnalysisManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
    }            
    finalResult = finalResult && result;
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()
        ->Message("delete", "empty file", fFileManager->GetFullFileName());
#endif
  }
  else {
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()
        ->Message("close", "file", fFileManager->GetFullFileName());
#endif
  }

  return finalResult;
} 
