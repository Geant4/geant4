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

// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#include "G4XmlAnalysisReader.hh"
#include "G4XmlRFileManager.hh"
#include "G4XmlRNtupleManager.hh"
#include "G4AnalysisUtilities.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4XmlAnalysisReader* G4XmlAnalysisReader::Instance()
{
  static G4ThreadLocalSingleton<G4XmlAnalysisReader> instance;
  return instance.Instance();
}

//_____________________________________________________________________________
G4XmlAnalysisReader::G4XmlAnalysisReader()
 : G4ToolsAnalysisReader("Xml")
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;

  // Create managers
  fNtupleManager = std::make_shared<G4XmlRNtupleManager>(fState);
  fFileManager = std::make_shared<G4XmlRFileManager>(fState);
  fNtupleManager->SetFileManager(fFileManager);

  // Set managers to base class
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4XmlAnalysisReader::~G4XmlAnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
}

//
// private methods
//

//_____________________________________________________________________________
G4bool G4XmlAnalysisReader::Reset()
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
G4bool  G4XmlAnalysisReader::CloseFilesImpl(G4bool reset)
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
