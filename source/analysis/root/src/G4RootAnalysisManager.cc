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
// $Id: G4RootAnalysisManager.cc 81633 2014-06-04 09:00:33Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4RootAnalysisManager.hh"
#include "G4RootFileManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4RootNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "tools/wroot/to"
#include "tools/wroot/file"

#include <iostream>
#include <cstdio>

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH1Mutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when merging H1 histograms 
  G4Mutex mergeH2Mutex = G4MUTEX_INITIALIZER;
}  

G4RootAnalysisManager* G4RootAnalysisManager::fgMasterInstance = 0;
G4ThreadLocal G4RootAnalysisManager* G4RootAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4RootAnalysisManager* G4RootAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4RootAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4RootAnalysisManager::G4RootAnalysisManager(G4bool isMaster)
 : G4VAnalysisManager("Root", isMaster),
   fH1Manager(0),
   fH2Manager(0),
   fNtupleManager(0),
   fFileManager(0)
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
  fH1Manager = new G4H1ToolsManager(fState);
  fH2Manager = new G4H2ToolsManager(fState);
  fNtupleManager = new G4RootNtupleManager(fState);
  fFileManager = new G4RootFileManager(fState);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4RootAnalysisManager::~G4RootAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = 0;
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH1()
{
  const std::vector<tools::histo::h1d*>& h1Vector
    = fH1Manager->GetH1Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fH1Manager->GetHnVector();

  if ( ! h1Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {
  
    for ( G4int i=0; i<G4int(h1Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->fActivation;
      G4String name = info->fName;
      // skip writing if activation is enabled and H1 is inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue; 
      tools::histo::h1d* h1 = h1Vector[i];
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "h1d", name);
#endif
      tools::wroot::directory* histoDirectory
        = fFileManager->GetHistoDirectory(); 
      G4bool result
        = to(*histoDirectory, *h1, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " <<  name << " failed";
        G4Exception("G4RootAnalysisManager::Write()",
                  "Analysis_W003", JustWarning, description);
        return false;       
      } 
    }
  }
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH1(&mergeH1Mutex);
    fgMasterInstance->fH1Manager->AddH1Vector(h1Vector);
    lH1.unlock();
  }  
  
  return true;
}
    
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::WriteH2()
{
  const std::vector<tools::histo::h2d*>& h2Vector
    = fH2Manager->GetH2Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fH2Manager->GetHnVector();

  if ( ! h2Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {
  
    for ( G4int i=0; i<G4int(h2Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->fActivation;
      G4String name = info->fName;
      // skip writing if inactivated
      if ( fState.GetIsActivation() && ( ! activation  ) ) continue;
      tools::histo::h2d* h2 = h2Vector[i];
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "h2d", name);
#endif
      tools::wroot::directory* histoDirectory
        = fFileManager->GetHistoDirectory(); 
      G4bool result
        = to(*histoDirectory, *h2, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4RootAnalysisManager::Write()",
                  "Analysis_W003", JustWarning, description);
        return false;       
      } 
    }
  }
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH2(&mergeH2Mutex);
    fgMasterInstance->fH2Manager->AddH2Vector(h2Vector);
    lH2.unlock();
  }
  
  return true;  
}  

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;
  
  G4bool result = fH1Manager->Reset();
  finalResult = finalResult && result;

  result = fH2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;
  G4bool result = fFileManager->SetFileName(fileName);
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

  G4bool finalResult = true;

  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) ) ) {
    G4ExceptionDescription description;
    description 
      << "      " << "No master G4RootAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram data will not be merged.";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W014", JustWarning, description);
  }
  
  // H1
  G4bool result = WriteH1();
  finalResult = finalResult && result;

  // H2
  result = WriteH2();
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
  G4bool result = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("close", "file", fFileManager->GetFullFileName());
#endif

  // reset data
  result = Reset();
  if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W002", JustWarning, description);
  } 

  // close file
  fFileManager->CloseFile();  

  // No files clean-up in sequential mode
  // or if ntuples are in use
  // (to avoid removing non-empty files if a sequential application
  // is built against MT Geant4 installation_
  if ( ! fState.IsMT() || ! fNtupleManager->IsEmpty() ) {
 #ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()
        ->Message("close", "file", fFileManager->GetFullFileName());
#endif
    return result;
  }  
  
  // Delete files if empty in MT mode
  if ( ( fState.GetIsMaster() && 
         fH1Manager->IsEmpty() && fH2Manager->IsEmpty() ) || 
       ( ( ! fState.GetIsMaster() ) && fNtupleManager->IsEmpty() ) ) {
    std::remove(fFileManager->GetFullFileName());
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

  return result;
} 
  
