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

// Author: Ivana Hrivnacova, 05/09/2014 (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisReader.hh"
#include "G4CsvRFileManager.hh"
#include "G4CsvRNtupleManager.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4CsvAnalysisReader* G4CsvAnalysisReader::Instance()
{
  static G4ThreadLocalSingleton<G4CsvAnalysisReader> instance;
  return instance.Instance();
}

//_____________________________________________________________________________
G4CsvAnalysisReader::G4CsvAnalysisReader()
 : G4ToolsAnalysisReader("Csv")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // Create managers
  fNtupleManager = std::make_shared<G4CsvRNtupleManager>(fState);
  fFileManager = std::make_shared<G4CsvRFileManager>(fState);
  fNtupleManager->SetFileManager(fFileManager);

  // Set managers to base class
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4CsvAnalysisReader::~G4CsvAnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
}

//
// private methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisReader::Reset()
{
// Reset histograms and ntuple

  auto result = true;

  result &= G4ToolsAnalysisReader::Reset();
  result &= fNtupleManager->Reset();

  return result;
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool  G4CsvAnalysisReader::CloseFilesImpl(G4bool reset)
{
  Message(kVL4, "close", "files");

  auto result = true;

  if (reset) {
    result &= Reset();
  }

  fFileManager->CloseFiles();

  Message(kVL2, "close", "files", "", result);

  return result;
}
