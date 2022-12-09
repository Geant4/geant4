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

#include "G4Hdf5FileManager.hh"
#include "G4Hdf5HnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4AutoLock.hh"

#include "toolx/hdf5/h2file"

using namespace G4Analysis;
using namespace tools;

//using namespace G4Analysis;

namespace {
  //Mutex to lock master manager when closing a file
  G4Mutex closeFileMutex = G4MUTEX_INITIALIZER;
}

//_____________________________________________________________________________
G4Hdf5FileManager::G4Hdf5FileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager<G4Hdf5File>(state)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4Hdf5HnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4Hdf5HnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4Hdf5HnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4Hdf5HnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4Hdf5HnFileManager<histo::p2d>>(this);
}

//
// private methods
//

//_____________________________________________________________________________
hid_t G4Hdf5FileManager::CreateDirectory(hid_t& file,
  const G4String& directoryName, const G4String& objectType)
{
// Method for both histograms and ntuples directories.

  // return if no file provided
  if (file < 0) return kInvalidId;

  // use default directory name if not provided
  auto newDirectoryName = directoryName;
  if ( newDirectoryName == "" ) {
      newDirectoryName = fgkDefaultDirectoryName;
      newDirectoryName += "_";
      newDirectoryName += objectType;
  }

  Message(kVL4, "create", "directory for " + objectType, newDirectoryName);

  auto success = true;

  // create directory
  auto directory = toolx_H5Gcreate(file, newDirectoryName, 0);
       // 0 seems to be an optional parameter. The web doc does not say what should
       // be the default value but 0 is what is found in examples, and in the code, if we pass 0, clearly some
       // default value is taken.
  if ( directory < 0 ) {
    Warn("Cannot create directory " + directoryName,
      fkClass, "CreateDirectory");
    success = false;
  }
  else {
    // write atb (header?)
    auto result = toolx::hdf5::write_atb(directory, "type", "directory");
    if ( !result) {
      Warn("Write_atb class failed for " + directoryName,
        fkClass, "CreateDirectory");
      success = false;
    }
  }

  Message(kVL2, "create", "directory for " + objectType, newDirectoryName, success);

  return directory;
}

//_____________________________________________________________________________
G4String G4Hdf5FileManager::GetNtupleFileName(Hdf5NtupleDescription* ntupleDescription)
{
  // get ntuple file name
  auto ntupleFileName = ntupleDescription->GetFileName();
  if (ntupleFileName.size() != 0u) {
    // update filename per object per thread
    ntupleFileName = GetTnFileName(ntupleFileName, GetFileType());
  }
  else {
    // get default file name
    ntupleFileName = GetFullFileName();
  }
  return ntupleFileName;
}

//
// protected methods
//

//_____________________________________________________________________________
std::shared_ptr<G4Hdf5File> G4Hdf5FileManager::CreateFileImpl(const G4String& fileName)
{
  // create a new file
  hid_t file = ::H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Do nothing if there is no file
  // (the error should be handled by caller)
  if ( file < 0 ) {
    Warn("::H5Fcreate failed " + fileName, fkClass, "CreateFileImpl");
    return std::make_shared<G4Hdf5File>(-1, -1, -1);
  }

  // create a header with general infos
  if(!toolx::hdf5::write_header(file)) {
    Warn("toolx::hdf5::write_header() failed for " + fileName,
      fkClass, "CreateFileImpl");
    return std::make_shared<G4Hdf5File>(-1, -1, -1);
  }

  // create histo directory
  auto hdirectory
    = CreateDirectory(file, fHistoDirectoryName, "histograms");
  if ( hdirectory < 0 ) {
    // Warning is issued in CreateDirectory
    return std::make_shared<G4Hdf5File>(-1, -1, -1);
  }

  // create ntuple directory
  auto ndirectory
    = CreateDirectory(file, fNtupleDirectoryName, "ntuples");
  if ( ndirectory < 0 ) {
    // Warnin is issued in CreateDirectory
    return std::make_shared<G4Hdf5File>(-1, -1, -1);
  }

  return std::make_shared<G4Hdf5File>(file, hdirectory, ndirectory);
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::WriteFileImpl(std::shared_ptr<G4Hdf5File> /*file*/)
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::CloseFileImpl(std::shared_ptr<G4Hdf5File> file)
{
  if ( ! file ) return false;

  G4AutoLock lock(&closeFileMutex);

  ::H5Gclose(std::get<1>(*file));
  ::H5Gclose(std::get<2>(*file));
  ::H5Fclose(std::get<0>(*file));

  lock.unlock();

  return true;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName = fileName;
  auto name = GetFullFileName();

  if ( fFile ) {
    Warn("File " + fileName + " already exists.", fkClass, "OpenFile");
    fFile.reset();
  }

  // create new file
  fFile = CreateTFile(name);
  if ( ! fFile ) {
    Warn("Failed to create file " + fileName, fkClass, "OpenFile");
    return false;
  }

  LockDirectoryNames();
  fIsOpenFile = true;

  return true;
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::CreateNtupleFile(
  Hdf5NtupleDescription* ntupleDescription)
{
  // get ntuple file name per object
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

  auto file = GetTFile(ntupleFileName, false);
  if (! file) {
    file = CreateTFile(ntupleFileName);
  }
  ntupleDescription->SetFile(file);

  return (ntupleDescription->GetFile() != nullptr);
}

//_____________________________________________________________________________
G4bool G4Hdf5FileManager::CloseNtupleFile(
  Hdf5NtupleDescription* ntupleDescription)
{
  // Notify not empty file
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);
  auto result = SetIsEmpty(ntupleFileName, ! ntupleDescription->GetHasFill());

  // Ntuple files are registered in file manager map.
  // they will be closed with CloseFiles() calls
  ntupleDescription->GetFile().reset();

  return result;
}

//_____________________________________________________________________________
hid_t G4Hdf5FileManager::GetHistoDirectory() const
{
  if ( ! fFile ) return kInvalidId;

  return std::get<1>(*fFile);
}

//_____________________________________________________________________________
hid_t G4Hdf5FileManager::GetNtupleDirectory() const
{
  if ( ! fFile ) return kInvalidId;

  return std::get<2>(*fFile);
}
