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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4RootFileManager.hh"
#include "G4RootHnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/wroot/file"
#include "tools/wroot/to"
#include "toolx/zlib"

using namespace tools;
using namespace G4Analysis;

//_____________________________________________________________________________
G4RootFileManager::G4RootFileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager<G4RootFile>(state)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4RootHnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4RootHnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4RootHnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4RootHnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4RootHnFileManager<histo::p2d>>(this);
}

//
// private methods
//

//_____________________________________________________________________________
tools::wroot::directory*  G4RootFileManager::CreateDirectory(
  tools::wroot::file* rfile,
  const G4String& directoryName,  [[maybe_unused]] const G4String& objectType) const
{
  if (rfile == nullptr) return nullptr;

  if ( directoryName == "" ) {
    // Do not create a new directory if its name is not set
    return &(rfile->dir());
  }

  Message(kVL4, "create", "directory for " + objectType, directoryName);

  auto directory = rfile->dir().mkdir(directoryName);
  if (directory == nullptr) {
    Warn("Cannot create directory " + directoryName, fkClass, "CreateDirectory");
    return nullptr;
  }
  Message(kVL2, "create", "directory for " + objectType, directoryName);

  return directory;
}

//_____________________________________________________________________________
G4String G4RootFileManager::GetNtupleFileName(
                              RootNtupleDescription* ntupleDescription,
                              G4bool perThread,
                              G4int mainNumber) const
{
  // get ntuple file name

  auto ntupleFileName = ntupleDescription->GetFileName();
  if (ntupleFileName.size() != 0u) {
    if ( perThread ) {
      ntupleFileName = GetTnFileName(ntupleFileName, GetFileType());
    }
  }
  else {
    // get default file name
    ntupleFileName = GetFullFileName(fFileName, perThread);
  }

  // update filename per mainNumber
  if ( mainNumber > -1 ) {
    // update filename per mainNumber
    ntupleFileName
      = G4Analysis::GetNtupleFileName(ntupleFileName, GetFileType(), mainNumber);
  }

  return ntupleFileName;
}

//
// protected methods
//

//_____________________________________________________________________________
std::shared_ptr<G4RootFile>
G4RootFileManager::CreateFileImpl(const G4String& fileName)
{
  // create file
  std::shared_ptr<wroot::file> file = std::make_shared<wroot::file>(G4cout, fileName);
  file->add_ziper('Z',toolx::compress_buffer);
  file->set_compression(fState.GetCompressionLevel());

  if ( ! file->is_open() ) {
    Warn("Cannot create file " + fileName, fkClass, "CreateFileImpl");
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  }

  // create histo directory
  tools::wroot::directory* hdirectory
    = CreateDirectory(file.get(), fHistoDirectoryName, "histograms");
  if (hdirectory == nullptr) {
    // Warning is issued in CreateDirectory
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  }

  // create ntuple directory
  tools::wroot::directory* ndirectory
    = CreateDirectory(file.get(), fNtupleDirectoryName, "ntuples");
  if (ndirectory == nullptr) {
    // Warning is issued in CreateDirectory
    return std::make_shared<G4RootFile>(nullptr, nullptr, nullptr);
  }

  return std::make_shared<G4RootFile>(file, hdirectory, ndirectory);
}

//_____________________________________________________________________________
G4bool G4RootFileManager::WriteFileImpl(std::shared_ptr<G4RootFile> file)
{
// New prototype: called by G4TFileManager base classe

  // nothing to be done, but file should exist
  return (file == nullptr) ? false : true;

}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseFileImpl(std::shared_ptr<G4RootFile> file)
{
// New prototype: called by G4TFileManager base classe

  if (file == nullptr) return false;

  // write file (only once)
  unsigned int n;
  std::get<0>(*file)->write(n);

  // close file
  std::get<0>(*file)->close();

  return true;
}

//
// public methods
//
//_____________________________________________________________________________
G4bool G4RootFileManager::OpenFile(const G4String& fileName)
{
// Open default file

  // Keep file name
  fFileName =  fileName;
  auto name = GetFullFileName();

  if ( fFile ) {
    Warn("File " + fileName + " already exists.", fkClass, "OpenFile");
    fFile.reset();
  }

  // Create file (and save in in the file map)
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
std::shared_ptr<G4RootFile> G4RootFileManager::CreateNtupleFile(
  RootNtupleDescription* ntupleDescription, G4int mainNumber)
{
  // get ntuple file name per object
  auto perThread = true;
  auto ntupleFileName = GetNtupleFileName(ntupleDescription, perThread, mainNumber);

  auto file = GetTFile(ntupleFileName, false);
  if (! file) {
    file = CreateTFile(ntupleFileName);
  }

  // register file in ntuple description only if it is not main ntuple file
  // (main ntuple files are managed in main ntuple manager)
  if ( mainNumber == -1 ) {
    ntupleDescription->SetFile(file);
  }

  return file;
}

//_____________________________________________________________________________
std::shared_ptr<G4RootFile> G4RootFileManager::GetNtupleFile(
  RootNtupleDescription* ntupleDescription,  G4bool perThread, G4int mainNumber) const
{
  // get ntuple file name per object
  auto ntupleFileName = GetNtupleFileName(ntupleDescription, perThread, mainNumber);

  return GetTFile(ntupleFileName, false);
}

//_____________________________________________________________________________
G4bool G4RootFileManager::CloseNtupleFile(
  RootNtupleDescription* ntupleDescription,  G4int mainNumber)
{
  // auto result = true;

  // Notify not empty file
  auto ntupleFileName = GetNtupleFileName(ntupleDescription, true, mainNumber);
  auto result = SetIsEmpty(ntupleFileName, ! ntupleDescription->GetHasFill());

  // Ntuple files will be closed with CloseFiles() calls
  ntupleDescription->GetFile().reset();

  return result;
}
