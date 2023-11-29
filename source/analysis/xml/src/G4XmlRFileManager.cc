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

#include "G4XmlRFileManager.hh"
#include "G4XmlHnRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "toolx/raxml"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4XmlRFileManager::G4XmlRFileManager(const G4AnalysisManagerState& state)
 : G4VRFileManager(state)
{
  // Create helpers defined in the base class
  fH1RFileManager = std::make_shared<G4XmlHnRFileManager<histo::h1d>>(this);
  fH2RFileManager = std::make_shared<G4XmlHnRFileManager<histo::h2d>>(this);
  fH3RFileManager = std::make_shared<G4XmlHnRFileManager<histo::h3d>>(this);
  fP1RFileManager = std::make_shared<G4XmlHnRFileManager<histo::p1d>>(this);
  fP2RFileManager = std::make_shared<G4XmlHnRFileManager<histo::p2d>>(this);
}

//_____________________________________________________________________________
G4XmlRFileManager::~G4XmlRFileManager()
{
  for ( auto& rfile : fRFiles ) {
    delete rfile.second;
  }

  delete fReadFactory;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlRFileManager::OpenRFile(const G4String& fileName)
{
  // Get full file name (add only extension)
  G4bool isPerThread = false;
  G4String name = GetFullFileName(fileName, isPerThread);

  Message(kVL4, "open", "read analysis file", name);

  G4bool verbose = false;

  // create factory (if it does not yet exist)
  if (fReadFactory == nullptr) {
    fReadFactory = new tools::xml::default_factory();
  }

  // create new file
  auto newFile = new toolx::raxml(*fReadFactory, G4cout, verbose);

  // clear objects
  // (this should not be needed when starting a new raxml)
  std::vector<tools::raxml_out>& objs = newFile->objects();
  objs.clear();

  G4bool compressed = false;
  if ( ! newFile->load_file(name, compressed) ) {
    Warn(G4String( "Cannot open file ") + name, fkClass, "OpenRFile");
    delete newFile;
    return false;
  }

  // add file in a map and delete the previous file if it exists
  auto it = fRFiles.find(name);
  if ( it != fRFiles.end() ) {
    delete it->second;
    it->second = newFile;
  }
  else {
    fRFiles[name] = newFile;
  }

  Message(kVL1, "open", "read analysis file", name);

  return true;
}

//_____________________________________________________________________________
toolx::raxml* G4XmlRFileManager::GetRFile(const G4String& fileName) const
{
  // Get full file name (add only extension)
  G4bool isPerThread = false;
  G4String name = GetFullFileName(fileName, isPerThread);

  auto it = fRFiles.find(name);
  if (it != fRFiles.end()) {
    return it->second;
  }
  return nullptr;
}
