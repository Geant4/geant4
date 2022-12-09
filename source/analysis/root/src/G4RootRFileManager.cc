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

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#include "G4RootRFileManager.hh"
#include "G4RootRFileDef.hh"
#include "G4RootHnRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/rroot/file"
#include "toolx/zlib"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4RootRFileManager::G4RootRFileManager(const G4AnalysisManagerState& state)
 : G4VRFileManager(state)
{
  // Create helpers defined in the base class
  fH1RFileManager = std::make_shared<G4RootHnRFileManager<histo::h1d>>(this);
  fH2RFileManager = std::make_shared<G4RootHnRFileManager<histo::h2d>>(this);
  fH3RFileManager = std::make_shared<G4RootHnRFileManager<histo::h3d>>(this);
  fP1RFileManager = std::make_shared<G4RootHnRFileManager<histo::p1d>>(this);
  fP2RFileManager = std::make_shared<G4RootHnRFileManager<histo::p2d>>(this);
}

//_____________________________________________________________________________
G4RootRFileManager::~G4RootRFileManager()
{
  // Delete all open file and their directories
  for ( auto& mapElement : fRFiles ) {
    auto rfileTuple = mapElement.second;
    delete std::get<1>(*rfileTuple);  // histo directory
    delete std::get<2>(*rfileTuple);; // ntuple directory
    delete std::get<0>(*rfileTuple);  // rfile
    delete rfileTuple;
  }
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4RootRFileManager::OpenRFile(const G4String& fileName,
                                     G4bool isPerThread)
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  Message(kVL4, "open", "read analysis file", name);

  // create new file
  auto newFile = new tools::rroot::file(G4cout, name);
  newFile->add_unziper('Z',toolx::decompress_buffer);

  if ( ! newFile->is_open() ) {
    Warn("Cannot open file " + name, fkClass, "OpenRFile");
    delete newFile;
    return false;
  }

  auto newFileTuple = new G4RootRFile(newFile, nullptr, nullptr);

  // add file in a map and delete the previous file if it exists
  auto it = fRFiles.find(name);
  if ( it != fRFiles.end() ) {
    delete it->second;
    it->second = newFileTuple;
  }
  else {
    fRFiles[name] = newFileTuple;
  }

  Message(kVL1, "open", "read analysis file", name);

  return true;
}

//_____________________________________________________________________________
G4RootRFile* G4RootRFileManager::GetRFile(const G4String& fileName,
                                          G4bool isPerThread) const
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  auto it = fRFiles.find(name);
  if (it != fRFiles.end()) {
    return it->second;
  }
  return nullptr;
}
