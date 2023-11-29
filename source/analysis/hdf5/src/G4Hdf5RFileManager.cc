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

#include "G4Hdf5RFileManager.hh"
#include "G4Hdf5HnRFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "toolx/hdf5/h2file"
#include "toolx/hdf5/group_exists"
#include "toolx/zlib"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4Hdf5RFileManager::G4Hdf5RFileManager(const G4AnalysisManagerState& state)
 : G4VRFileManager(state)
{
  // Create helpers defined in the base class
  fH1RFileManager = std::make_shared<G4Hdf5HnRFileManager<histo::h1d>>(this);
  fH2RFileManager = std::make_shared<G4Hdf5HnRFileManager<histo::h2d>>(this);
  fH3RFileManager = std::make_shared<G4Hdf5HnRFileManager<histo::h3d>>(this);
  fP1RFileManager = std::make_shared<G4Hdf5HnRFileManager<histo::p1d>>(this);
  fP2RFileManager = std::make_shared<G4Hdf5HnRFileManager<histo::p2d>>(this);
}

//
// private methods
//

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::OpenRFile(const G4String& fileName,
                                     G4bool isPerThread)
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  Message(kVL4, "open", "read analysis file", name);

  // create new file
  hid_t newFile = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( newFile < 0 ) {
    Warn("Cannot open file " + name, fkClass, "OpenRFile");
    return kInvalidId;
  }

  // newFile->add_unziper('Z',toolx::decompress_buffer);

  // add file in a map
  fRFiles[name] = G4Hdf5File(newFile, kInvalidId, kInvalidId);

  Message(kVL1, "open", "read analysis file", name);

  return newFile;
}

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::OpenDirectory(hid_t file, const G4String& directoryName)
{
  Message(kVL4, "open", "read directory", directoryName);

  auto directory = toolx_H5Gopen(file, directoryName);
  if ( directory < 0 ) {
    Warn("Cannot open directory " + directoryName, fkClass, "OpenDirectory");
    return kInvalidId;
  }
  Message(kVL2, "open", "read directory", directoryName);
  return directory;
}

//_____________________________________________________________________________
hid_t  G4Hdf5RFileManager::GetRDirectory(const G4String& directoryType,
                                         const G4String& fileName,
                                         const G4String& dirName,
                                         G4bool isPerThread)
{
  // Get or open a file
  auto rfile = GetRFile(fileName, isPerThread);
  if (rfile == nullptr) {
    // Try to open it if not found in the map
    if ( OpenRFile(fileName, isPerThread) < 0 ) return kInvalidId;
    rfile = GetRFile(fileName, isPerThread);
  }

  auto isHistograms = (directoryType == "histograms");

  // Get directory if already open
  hid_t directory = kInvalidId;
  if ( isHistograms ) {
    directory = std::get<1>(*rfile);
  } else {
    directory = std::get<2>(*rfile);
  }
  if ( directory != kInvalidId ) {
    return directory;
  }

  // Use default directory name if not specified
  auto newDirName = dirName;
  if ( newDirName == "" ) {
      // Create the default directory if the name is not set and the default directory
      // does not yet exist
      newDirName = fgkDefaultDirectoryName;
      newDirName += "_";
      newDirName += directoryType;
  }

  // Open directory
  directory = OpenDirectory(std::get<0>(*rfile), newDirName);

  // Update
  if ( isHistograms ) {
    std::get<1>(*rfile) = directory;
  } else {
    std::get<2>(*rfile) = directory;
  }

  return directory;
}

//
// public methods
//

//_____________________________________________________________________________
G4Hdf5File* G4Hdf5RFileManager::GetRFile(const G4String& fileName,
                                         G4bool isPerThread)
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  auto it = fRFiles.find(name);
  if (it != fRFiles.end()) {
    return &(it->second);
  }
  return nullptr;
}

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::GetHistoRDirectory(const G4String& fileName, const G4String& dirName,
                                             G4bool isPerThread)
{
  return GetRDirectory("histograms", fileName, dirName, isPerThread);
}

//_____________________________________________________________________________
hid_t G4Hdf5RFileManager::GetNtupleRDirectory(const G4String& fileName, const G4String& dirName,
                                              G4bool isPerThread)
{
  return GetRDirectory("ntuples", fileName, dirName, isPerThread);
}

//_____________________________________________________________________________
void G4Hdf5RFileManager::CloseFiles()
{
  // Close all open directories and file
  for ( auto [key, rfile] : fRFiles ) {
    if (std::get<1>(rfile) != kInvalidId) {
      ::H5Gclose(std::get<1>(rfile));
    }
    if (std::get<2>(rfile) != kInvalidId) {
      ::H5Gclose(std::get<2>(rfile));
    }
    if (std::get<0>(rfile) != kInvalidId) {
      ::H5Fclose(std::get<0>(rfile));
    }
  }
}
