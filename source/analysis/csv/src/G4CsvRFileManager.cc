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

// Author: Ivana Hrivnacova, 21/10/2014  (ivana@ipno.in2p3.fr)

#include "G4CsvRFileManager.hh"
#include "G4CsvHnRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4CsvRFileManager::G4CsvRFileManager(const G4AnalysisManagerState& state)
 : G4VRFileManager(state)
{
  // Create helpers defined in the base class
  fH1RFileManager = std::make_shared<G4CsvHnRFileManager<histo::h1d>>(this);
  fH2RFileManager = std::make_shared<G4CsvHnRFileManager<histo::h2d>>(this);
  fH3RFileManager = std::make_shared<G4CsvHnRFileManager<histo::h3d>>(this);
  fP1RFileManager = std::make_shared<G4CsvHnRFileManager<histo::p1d>>(this);
  fP2RFileManager = std::make_shared<G4CsvHnRFileManager<histo::p2d>>(this);
}

//_____________________________________________________________________________
G4CsvRFileManager::~G4CsvRFileManager()
{
  for ( auto& rfile : fRFiles ) {
    delete rfile.second;
  }
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvRFileManager::OpenRFile(const G4String& fileName)
{
  Message(kVL4, "open", "read analysis file", fileName);

  // create new file
  auto newFile = new std::ifstream(fileName);
  if ( ! newFile->is_open() ) {
    Warn("Cannot open file " + fileName, fkClass, "OpenRFile");
    return false;
  }

  // add file in a map and delete the previous file if it exists
  auto it = fRFiles.find(fileName);
  if ( it != fRFiles.end() ) {
    delete it->second;
    it->second = newFile;
  }
  else {
    fRFiles[fileName] = newFile;
  }

  Message(kVL1, "open", "read analysis file", fileName);

  return true;
}

//_____________________________________________________________________________
std::ifstream* G4CsvRFileManager::GetRFile(const G4String& fileName) const
{
  auto it = fRFiles.find(fileName);
  if (it != fRFiles.end()) {
    return it->second;
  }
  return nullptr;
}
