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
#include "G4CsvNtupleFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4UnitsTable.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4CsvAnalysisManager* G4CsvAnalysisManager::Instance()
{
  static G4ThreadLocalSingleton<G4CsvAnalysisManager> instance;
  fgIsInstance = true;
  return instance.Instance();
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::IsInstance()
{
  return fgIsInstance;
}

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager()
 : G4ToolsAnalysisManager("Csv")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // File Manager
  fFileManager = std::make_shared<G4CsvFileManager>(fState);
  SetFileManager(fFileManager);

  // Ntuple file manager
  fNtupleFileManager = std::make_shared<G4CsvNtupleFileManager>(fState);
  fNtupleFileManager->SetFileManager(fFileManager);
  fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgIsInstance = false;
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFileImpl(const G4String& fileName)
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
G4bool G4CsvAnalysisManager::WriteImpl()
{
  // nothing to be done for Csv file
  auto result = true;

  Message(kVL4, "write", "files");

  if ( G4Threading::IsWorkerThread() )  {
    result &= G4ToolsAnalysisManager::Merge();
  }
  else {
    // Write all histograms/profile on master
    result &= G4ToolsAnalysisManager::WriteImpl();
  }

  // Ntuples
  // Nothing to be done

  // Write ASCII if activated
  // Not available
  //if ( IsAscii() ) {
  //  result &= WriteAscii();
  //}

  Message(kVL3, "write", "files", "", result);

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFileImpl(G4bool reset)
{
  auto result = true;

  // close open files
  result &= fFileManager->CloseFiles();

  // Histogram and profile files are created/closed indivudually for each
  // object in WriteT function

  result &= fNtupleFileManager->ActionAtCloseFile(reset);

  // Reset data
  if ( reset ) {
    result = Reset();
    if ( ! result ) {
      Warn("Resetting data failed", fkClass, "CloseFileImpl");
    }
  }

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::ResetImpl()
{
// Reset histograms and ntuple

  auto result = true;

  result &= G4ToolsAnalysisManager::ResetImpl();
  if ( fNtupleFileManager != nullptr ) {
    result &= fNtupleFileManager->Reset();
  }

  return result;
}
