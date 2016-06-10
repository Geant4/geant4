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
// $Id: G4CsvAnalysisManager.cc 85317 2014-10-27 15:58:15Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisManager.hh"
#include "G4CsvFileManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4CsvNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisManagerState.hh"
#include "G4UnitsTable.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "tools/wcsv_histo"

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
   fH1Manager(0),
   fH2Manager(0),
   fH3Manager(0),
   fP1Manager(0),
   fP2Manager(0),
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
  fH1Manager = new G4H1ToolsManager(fState);
  fH2Manager = new G4H2ToolsManager(fState);
  fH3Manager = new G4H3ToolsManager(fState);
  fP1Manager = new G4P1ToolsManager(fState);
  fP2Manager = new G4P2ToolsManager(fState);
  fNtupleManager = new G4CsvNtupleManager(fState);
  fFileManager = new G4CsvFileManager(fState);
  fNtupleManager->SetFileManager(fFileManager);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetH1Manager(fH1Manager);
  SetH2Manager(fH2Manager);
  SetH3Manager(fH3Manager);
  SetP1Manager(fP1Manager);
  SetP2Manager(fP2Manager);
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
    = fNtupleManager->GetNtupleDescriptionVector();

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

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteH1()
{
  const std::vector<tools::histo::h1d*>& h1Vector
    = fH1Manager->GetH1Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fH1Manager->GetHnVector();

  if ( ! h1Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {

    for ( G4int i=0; i<G4int(h1Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->GetActivation();
      G4String name = info->GetName();
      // skip writing if activation is enabled and H1 is inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue; 
      tools::histo::h1d* h1 = h1Vector[i];
      G4String fileName = fFileManager->GetHnFileName("h1", name);
      std::ofstream hnFile(fileName);

      G4bool result
        = tools::wcsv::hto(hnFile, h1->s_cls(), *h1);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4CsvAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
      hnFile.close();
#ifdef G4VERBOSE
      if ( fState.GetVerboseL1() ) 
        fState.GetVerboseL1()->Message("write", "file", fileName);
#endif
      fFileManager->LockHistoDirectoryName();
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
G4bool G4CsvAnalysisManager::WriteH2()
{
  const std::vector<tools::histo::h2d*>& h2Vector
    = fH2Manager->GetH2Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fH2Manager->GetHnVector();

  if ( ! h2Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {

    // h2 histograms
    for ( G4int i=0; i<G4int(h2Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->GetActivation();
      G4String name = info->GetName();
      // skip writing if inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue;
      tools::histo::h2d* h2 = h2Vector[i];
      G4String fileName = fFileManager->GetHnFileName("h2", name);
      std::ofstream hnFile(fileName);
      G4bool result
        = tools::wcsv::hto(hnFile, h2->s_cls(), *h2);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4CsvAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
      hnFile.close();
#ifdef G4VERBOSE
      if ( fState.GetVerboseL1() ) 
        fState.GetVerboseL1()->Message("write", "file", fileName);
#endif
      fFileManager->LockHistoDirectoryName();
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
G4bool G4CsvAnalysisManager::WriteH3()
{
  const std::vector<tools::histo::h3d*>& h3Vector
    = fH3Manager->GetH3Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fH3Manager->GetHnVector();

  if ( ! h3Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {

    // h3 histograms
    for ( G4int i=0; i<G4int(h3Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->GetActivation();
      G4String name = info->GetName();
      // skip writing if inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue;
      tools::histo::h3d* h3 = h3Vector[i];
      G4String fileName = fFileManager->GetHnFileName("h3", name);
      std::ofstream hnFile(fileName);
      G4bool result
        = tools::wcsv::hto(hnFile, h3->s_cls(), *h3);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4CsvAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
      hnFile.close();
#ifdef G4VERBOSE
      if ( fState.GetVerboseL1() ) 
        fState.GetVerboseL1()->Message("write", "file", fileName);
#endif
      fFileManager->LockHistoDirectoryName();
    }
  }  
  else {
    // The worker manager just adds its histograms to the master
    // This operation needs a lock
    G4AutoLock lH3(&mergeH3Mutex);
    fgMasterInstance->fH3Manager->AddH3Vector(h3Vector);
    lH3.unlock();
  }  
  
  return true;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteP1()
{
  const std::vector<tools::histo::p1d*>& p1Vector
    = fP1Manager->GetP1Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fP1Manager->GetHnVector();

  if ( ! p1Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {
  
    for ( G4int i=0; i<G4int(p1Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->GetActivation();
      G4String name = info->GetName();
      // skip writing if activation is enabled and P1 is inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue; 
      tools::histo::p1d* p1 = p1Vector[i];
      G4String fileName = fFileManager->GetHnFileName("p1", name);
      std::ofstream hnFile(fileName);
      G4bool result
        = tools::wcsv::pto(hnFile, p1->s_cls(), *p1);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving profile " << name << " failed";
        G4Exception("G4CsvAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
      hnFile.close();
#ifdef G4VERBOSE
      if ( fState.GetVerboseL1() ) 
        fState.GetVerboseL1()->Message("write", "file", fileName);
#endif
      fFileManager->LockProfileDirectoryName();
    }
  }
  else {
    // The worker manager just adds its profiles to the master
    // This operation needs a lock
    G4AutoLock lP1(&mergeP1Mutex);
    fgMasterInstance->fP1Manager->AddP1Vector(p1Vector);
    lP1.unlock();
  }  
  
  return true;
}
    
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteP2()
{
  const std::vector<tools::histo::p2d*>& p2Vector
    = fP2Manager->GetP2Vector();
  const std::vector<G4HnInformation*>& hnVector
    = fP2Manager->GetHnVector();

  if ( ! p2Vector.size() ) return true;

  if ( ! G4Threading::IsWorkerThread() )  {
  
    for ( G4int i=0; i<G4int(p2Vector.size()); ++i ) {
      G4HnInformation* info = hnVector[i];
      G4bool activation = info->GetActivation();
      G4String name = info->GetName();
      // skip writing if activation is enabled and P2 is inactivated
      if ( fState.GetIsActivation() && ( ! activation ) ) continue; 
      tools::histo::p2d* p2 = p2Vector[i];
      G4String fileName = fFileManager->GetHnFileName("p2", name);
      std::ofstream hnFile(fileName);
      G4bool result
        = tools::wcsv::pto(hnFile, p2->s_cls(), *p2);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving profile " << name << " failed";
        G4Exception("G4CsvAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
      hnFile.close();
#ifdef G4VERBOSE
      if ( fState.GetVerboseL1() ) 
        fState.GetVerboseL1()->Message("write", "file", fileName);
#endif
      fFileManager->LockProfileDirectoryName();
    }
  }
  else {
    // The worker manager just adds its profiles to the master
    // This operation needs a lock
    G4AutoLock lP2(&mergeP2Mutex);
    fgMasterInstance->fP2Manager->AddP2Vector(p2Vector);
    lP2.unlock();
  }  
  
  return true;
}
    
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;

  G4bool result = fH1Manager->Reset();
  finalResult = finalResult && result;

  result = fH2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fH3Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP1Manager->Reset();
  finalResult = finalResult && result;
  
  result = fP2Manager->Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;
  G4bool result = fFileManager->SetFileName(fileName);
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
  G4bool finalResult = true;
  
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
  G4bool result = WriteH1();
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
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("write", "files", "", finalResult);
#endif

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFileImpl()
{
  G4bool finalResult = true;

  // Unlock file name only
  G4bool result = fFileManager->CloseFile();
  finalResult = finalResult && result;
   
  // Histogram and profile files are created/closed indivudually for each
  // object in WriteHn{Pn] function
 
  // Close ntuple files
  if ( ( ! G4Threading::IsMultithreadedApplication() ) || 
       ( ! fState.GetIsMaster() ) ) {
    // In sequential mode or in MT mode only on workers
    result = CloseNtupleFiles();
    finalResult = finalResult && result;
  }  

  // reset data
  result = Reset();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "Resetting data failed";
    G4Exception("G4CsvAnalysisManager::CloseFile()",
              "Analysis_W021", JustWarning, description);
    result = false;       
  } 
  finalResult = finalResult && result;

  return finalResult; 
} 
