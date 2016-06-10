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
// $Id: G4XmlAnalysisManager.cc 85317 2014-10-27 15:58:15Z gcosmo $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4XmlAnalysisManager.hh"
#include "G4XmlFileManager.hh"
#include "G4H1ToolsManager.hh"
#include "G4H2ToolsManager.hh"
#include "G4H3ToolsManager.hh"
#include "G4P1ToolsManager.hh"
#include "G4P2ToolsManager.hh"
#include "G4XmlNtupleManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "tools/waxml/histos"

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

G4XmlAnalysisManager* G4XmlAnalysisManager::fgMasterInstance = 0;
G4ThreadLocal G4XmlAnalysisManager* G4XmlAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4XmlAnalysisManager* G4XmlAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4XmlAnalysisManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4XmlAnalysisManager::G4XmlAnalysisManager(G4bool isMaster)
 : G4VAnalysisManager("Xml", isMaster),
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
    description 
      << "      " 
      << "G4XmlAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4XmlAnalysisManager::G4XmlAnalysisManager",
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
  fNtupleManager = new G4XmlNtupleManager(fState);
  fFileManager = new G4XmlFileManager(fState);
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
G4XmlAnalysisManager::~G4XmlAnalysisManager()
{  
  if ( fState.GetIsMaster() ) fgMasterInstance = 0;
  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::WriteH1()
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
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "h1d", name);
#endif
      G4String path = "/";
      path.append(fFileManager->GetHistoDirectoryName());
      std::ofstream* hnFile = fFileManager->GetHnFile();
      G4bool result
        = tools::waxml::write(*hnFile, *h1, path, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
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
G4bool G4XmlAnalysisManager::WriteH2()
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
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "h2d", name);
#endif
      G4String path = "/";
      path.append(fFileManager->GetHistoDirectoryName());
      std::ofstream* hnFile = fFileManager->GetHnFile();
      G4bool result
        = tools::waxml::write(*hnFile, *h2, path, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
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
G4bool G4XmlAnalysisManager::WriteH3()
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
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "h3d", name);
#endif
      G4String path = "/";
      path.append(fFileManager->GetHistoDirectoryName());
      std::ofstream* hnFile = fFileManager->GetHnFile();
      G4bool result
        = tools::waxml::write(*hnFile, *h3, path, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving histogram " << name << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
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
G4bool G4XmlAnalysisManager::WriteP1()
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
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "p1d", name);
#endif
      G4String path = "/";
      path.append(fFileManager->GetProfileDirectoryName());
      std::ofstream* hnFile = fFileManager->GetHnFile();
      G4bool result
        = tools::waxml::write(*hnFile, *p1, path, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving profile " << name << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
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
G4bool G4XmlAnalysisManager::WriteP2()
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
#ifdef G4VERBOSE
      if ( fState.GetVerboseL3() ) 
        fState.GetVerboseL3()->Message("write", "p2d", name);
#endif
      G4String path = "/";
      path.append(fFileManager->GetProfileDirectoryName());
      std::ofstream* hnFile = fFileManager->GetHnFile();
      G4bool result
        = tools::waxml::write(*hnFile, *p2, path, name);
      if ( ! result ) {
        G4ExceptionDescription description;
        description << "      " << "saving profile " << name << " failed";
        G4Exception("G4XmlAnalysisManager::Write()",
                  "Analysis_W022", JustWarning, description);
        return false;       
      } 
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
G4bool G4XmlAnalysisManager::WriteNtuple()
{
  const std::vector<G4XmlNtupleDescription*>& ntupleVector
    = fNtupleManager->GetNtupleDescriptionVector();

  for ( G4int i=0; i<G4int(ntupleVector.size()); ++i ) {
    if ( ntupleVector[i]->fNtuple ) ntupleVector[i]->fNtuple->write_trailer();
  }
  
  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseNtupleFiles()
{
  const std::vector<G4XmlNtupleDescription*>& ntupleVector
    = fNtupleManager->GetNtupleDescriptionVector();

  // Close ntuple files
  std::vector<G4XmlNtupleDescription*>::const_iterator it;  
  for (it = ntupleVector.begin(); it != ntupleVector.end(); it++ ) {
    fFileManager->CloseNtupleFile((*it));
  }
  
  return true;
}    


//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::Reset()
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
 
// 
// protected methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4bool finalResult = true;
  G4bool result = fFileManager->SetFileName(fileName);
  finalResult = finalResult && result;

#ifdef G4VERBOSE
  G4String name = fFileManager->GetFullFileName();
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("open", "analysis file", name);
  }  
#endif

  // Only lock file name in file manager
  result = fFileManager->OpenFile(fileName);
  finalResult = finalResult && result;

  // Create histograms file (on master)
  if ( fState.GetIsMaster() ) {
    result = fFileManager->CreateHnFile();
    finalResult = finalResult && result;
  }  

  // Create ntuples if they are booked
  // (The files will be created with creating ntuples)
  fNtupleManager->CreateNtuplesFromBooking();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("open", "analysis file", name);
#endif
  
  return finalResult;
}  
  
//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::WriteImpl() 
{
  G4bool finalResult = true;

#ifdef G4VERBOSE
  G4String name = fFileManager->GetFullFileName();
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write", "files", name);
#endif

  // ntuples 
  WriteNtuple();

  if ( ! fgMasterInstance && 
       ( ( ! fH1Manager->IsEmpty() ) || ( ! fH2Manager->IsEmpty() ) || 
         ( ! fH3Manager->IsEmpty() ) || ( ! fP1Manager->IsEmpty() ) ||
         ( ! fP2Manager->IsEmpty() ) ) ) {

    G4ExceptionDescription description;
    description 
      << "      " << "No master G4XmlAnalysisManager instance exists." 
      << G4endl 
      << "      " << "Histogram data will not be merged.";
      G4Exception("G4XmlAnalysisManager::Write()",
                "Analysis_W031", JustWarning, description);
                
    // Create Hn file per thread
    G4bool result = fFileManager->CreateHnFile();
    if ( ! result ) return false;       
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

  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii(fFileManager->GetFileName());
    finalResult = finalResult && result;
  }   

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("write", "file", fFileManager->GetFullFileName(), finalResult);
#endif

  return finalResult;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseFileImpl()
{
  G4bool finalResult = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "files", "");
#endif

  // Unlock file name only
  G4bool result = fFileManager->CloseFile();
  finalResult = finalResult && result;
  
  // Close Hn file
  result = fFileManager->CloseHnFile();  
  finalResult = finalResult && result;
  
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
      G4Exception("G4XmlAnalysisManager::CloseFile()",
                "Analysis_W021", JustWarning, description);
  } 
  finalResult = finalResult && result;

  // delete files if empty
  // (ntuple files are created only if an ntuple is created)
  if ( fFileManager->GetHnFile() && 
       fH1Manager->IsEmpty() && fH2Manager->IsEmpty() && fH3Manager->IsEmpty() &&
       fP1Manager->IsEmpty() && fP2Manager->IsEmpty() ) {
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
    if ( fState.GetVerboseL2() ) 
      fState.GetVerboseL2()
        ->Message("close", "files", "");
#endif
  }

  return finalResult; 
}
