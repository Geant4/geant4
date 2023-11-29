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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4Hdf5AnalysisManager.hh"
#include "G4Hdf5FileManager.hh"
#include "G4Hdf5NtupleFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

using namespace G4Analysis;

// mutex in a file scope

namespace {
  //Mutex to lock master manager when opening a file
  G4Mutex openFileMutex = G4MUTEX_INITIALIZER;
  //Mutex to lock master manager when closing a file
  G4Mutex closeFileMutex = G4MUTEX_INITIALIZER;
}

// G4Hdf5AnalysisManager* G4Hdf5AnalysisManager::fgMasterInstance = nullptr;
// G4ThreadLocal G4bool G4Hdf5AnalysisManager::fgIsInstance = false;

//_____________________________________________________________________________
G4Hdf5AnalysisManager* G4Hdf5AnalysisManager::Instance()
{
  static G4ThreadLocalSingleton<G4Hdf5AnalysisManager> instance;
  fgIsInstance = true;
  return instance.Instance();
}

//_____________________________________________________________________________
G4bool G4Hdf5AnalysisManager::IsInstance()
{
  return fgIsInstance;
}

//_____________________________________________________________________________
G4Hdf5AnalysisManager::G4Hdf5AnalysisManager()
 : G4ToolsAnalysisManager("Hdf5")
{
#ifdef G4MULTITHREADED
#ifndef H5_HAVE_THREADSAFE
    G4Exception("G4Hdf5AnalysisManager::G4Hdf5AnalysisManager",
                "Analysis_F001", FatalException,
                "Your HDF5 lib is not built with H5_HAVE_THREADSAFE.");
#endif
#endif

  // File manager
  auto fileManager = std::make_shared<G4Hdf5FileManager>(fState);
  SetFileManager(fileManager);

  // Ntuple file manager
  fNtupleFileManager = std::make_shared<G4Hdf5NtupleFileManager>(fState);
  SetNtupleFileManager(fNtupleFileManager);
  fNtupleFileManager->SetFileManager(fileManager);
  fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
}

//_____________________________________________________________________________
G4Hdf5AnalysisManager::~G4Hdf5AnalysisManager()
{
  fgIsInstance = false;
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4Hdf5AnalysisManager::OpenFileImpl(const G4String& fileName)
{
  G4AutoLock lock(&openFileMutex);
  auto result = G4ToolsAnalysisManager::OpenFileImpl(fileName);
  lock.unlock();

  return result;
}

//_____________________________________________________________________________
G4bool G4Hdf5AnalysisManager::CloseFileImpl(G4bool reset)
{
  G4AutoLock lock(&closeFileMutex);
  auto result = G4ToolsAnalysisManager::CloseFileImpl(reset);
  lock.unlock();

  return result;
}
