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
  // File Manager
  auto fileManager = std::make_shared<G4CsvFileManager>(fState);
  SetFileManager(fileManager);

  // Ntuple file manager
  fNtupleFileManager = std::make_shared<G4CsvNtupleFileManager>(fState);
  SetNtupleFileManager(fNtupleFileManager);
  fNtupleFileManager->SetFileManager(fileManager);
  fNtupleFileManager->SetBookingManager(fNtupleBookingManager);
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{
  fgIsInstance = false;
}
