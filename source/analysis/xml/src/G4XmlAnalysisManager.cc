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

#include "G4XmlAnalysisManager.hh"
#include "G4XmlFileManager.hh"
#include "G4XmlNtupleFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"

using namespace G4Analysis;
using std::make_shared;

//_____________________________________________________________________________
G4XmlAnalysisManager* G4XmlAnalysisManager::Instance()
{
  static G4ThreadLocalSingleton<G4XmlAnalysisManager> instance;
  fgIsInstance = true;
  return instance.Instance();
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::IsInstance()
{
  return fgIsInstance;
}

//_____________________________________________________________________________
G4XmlAnalysisManager::G4XmlAnalysisManager()
 : G4ToolsAnalysisManager("Xml")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // File Manager
  fFileManager = std::make_shared<G4XmlFileManager>(fState);
  SetFileManager(fFileManager);

  // Ntuple file manager
  fNtupleFileManager = std::make_shared<G4XmlNtupleFileManager>(fState);
  fNtupleFileManager->SetFileManager(fFileManager);
  fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
}

//_____________________________________________________________________________
G4XmlAnalysisManager::~G4XmlAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgIsInstance = false;
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::OpenFileImpl(const G4String& fileName)
{
  // Create ntuple manager(s)
  // and set it to base class which takes then their ownership
  SetNtupleManager(fNtupleFileManager->CreateNtupleManager());

  auto result = true;

  // Save file name in file manager
  result &= fFileManager->OpenFile(fileName);

  // Open ntuple files and create ntuples from bookings
  result &= fNtupleFileManager->ActionAtOpenFile(fFileManager->GetFullFileName());

  return result;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::WriteImpl()
{
  auto result = true;

  Message(kVL4, "write", "files");

  // ntuples
  fNtupleFileManager->ActionAtWrite();

  if ( G4Threading::IsWorkerThread() )  {
    result &= G4ToolsAnalysisManager::Merge();
  }
  else {
    // Open all files registered with objects
    fFileManager->OpenFiles();

    // Write all histograms/profile on master
    result &= G4ToolsAnalysisManager::WriteImpl();
  }

  // Write ASCII if activated
  if ( IsAscii() ) {
    result &= WriteAscii(fFileManager->GetFileName());
  }

  // File
  result &= fFileManager->WriteFiles();

  Message(kVL3, "write", "files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::CloseFileImpl(G4bool reset)
{
  auto result = true;

  Message(kVL4, "close", "files");

  // Close open files
  result &= fFileManager->CloseFiles();

  // Close ntuple files
  result &= fNtupleFileManager->ActionAtCloseFile(reset);

  // Reset data
  if ( reset ) {
    if ( ! Reset() ) {
      Warn("Resetting data failed", fkClass, "CloseFileImpl");
      result = false;
    }
  }

  // delete files if empty
  // (ntuple files are created only if an ntuple is created)
  if ( fFileManager->GetFile() && G4ToolsAnalysisManager::IsEmpty() ) {
    if ( std::remove(fFileManager->GetFullFileName()) ) {
      //  std::remove returns 0 when success
      Warn(G4String(
        "Removing file ") + fFileManager->GetFullFileName() + " failed",
        fkClass, "CloseFileImpl");
      result = false;
    }
    Message(kVL1, "delete", "empty file", fFileManager->GetFullFileName());
  }
  else {
    Message(kVL3, "close", "files");
  }

  return result;
}

//_____________________________________________________________________________
G4bool G4XmlAnalysisManager::ResetImpl()
{
// Reset histograms and ntuple

  auto result = true;

  result &= G4ToolsAnalysisManager::ResetImpl();
  if ( fNtupleFileManager != nullptr ) {
    result &= fNtupleFileManager->Reset();
  }

  return result;
}
