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

#include "G4Hdf5AnalysisReader.hh"
#include "G4Hdf5RFileManager.hh"
#include "G4Hdf5RNtupleManager.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4Hdf5AnalysisReader* G4Hdf5AnalysisReader::Instance()
{
  static G4ThreadLocalSingleton<G4Hdf5AnalysisReader> instance;
  return instance.Instance();
}

//_____________________________________________________________________________
G4Hdf5AnalysisReader::G4Hdf5AnalysisReader()
 : G4ToolsAnalysisReader("Hdf5")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // Create managers
  fNtupleManager = std::make_shared<G4Hdf5RNtupleManager>(fState);
  fFileManager = std::make_shared<G4Hdf5RFileManager>(fState);
  fNtupleManager->SetFileManager(fFileManager);

  // Set managers to base class
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4Hdf5AnalysisReader::~G4Hdf5AnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
}

//
// private methods
//

//_____________________________________________________________________________
G4bool G4Hdf5AnalysisReader::Reset()
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
G4bool  G4Hdf5AnalysisReader::CloseFilesImpl(G4bool reset)
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
