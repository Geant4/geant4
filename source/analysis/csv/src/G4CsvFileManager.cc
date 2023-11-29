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

#include "G4CsvFileManager.hh"
#include "G4CsvHnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Filesystem.hh"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4CsvFileManager::G4CsvFileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager(state)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4CsvHnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4CsvHnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4CsvHnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4CsvHnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4CsvHnFileManager<histo::p2d>>(this);
}

//
// private methods
//

//_____________________________________________________________________________
G4String G4CsvFileManager::GetNtupleFileName(CsvNtupleDescription* ntupleDescription)
{
  // get ntuple file name
  auto ntupleFileName = ntupleDescription->GetFileName();
  auto cycle = GetCycle();
  if (ntupleFileName.size() != 0u) {
    // update filename per object per thread
    ntupleFileName = GetTnFileName(ntupleFileName, GetFileType(), cycle);
  }
  else {
    // compose ntuple file name from the default file name
    ntupleFileName = GetNtupleFileName(ntupleDescription->GetNtupleBooking().name(), cycle);
  }

  if ( IsNtupleDirectory() ) {
    ntupleFileName = "./" + GetNtupleDirectoryName() + "/" + ntupleFileName;
  }

  return ntupleFileName;
}

//
// protected methods
//

//_____________________________________________________________________________
std::shared_ptr<std::ofstream> G4CsvFileManager::CreateFileImpl(const G4String& fileName)
{
  std::shared_ptr<std::ofstream> file = std::make_shared<std::ofstream>(fileName);
  if ( file->fail() ) {
    Warn("Cannot create file " + fileName, fkClass, "CreateFileImpl");
    return nullptr;
  }

  return file;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::WriteFileImpl(std::shared_ptr<std::ofstream> /*file*/)
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::CloseFileImpl(std::shared_ptr<std::ofstream> file)
{
  if ( ! file ) return false;

  // close file
  file->close();

  return true;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvFileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;

  fIsOpenFile = true;

  return true;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::SetHistoDirectoryName(const G4String& dirName)
{
  // A directory is taken into account only if it exists in file system
  if ( G4fs::is_directory(dirName.data()) ) {
     fIsHistoDirectory = G4VFileManager::SetHistoDirectoryName(dirName);
     return fIsHistoDirectory;
  }

  G4Analysis::Warn("Directory " + dirName + " does not exists.\n"
    "Histograms will be written in the current directory.",
    fkClass, "SetHistoDirectoryName");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::SetNtupleDirectoryName(const G4String& dirName)
{
  // A directory is taken into account only if it exists in file system
  if ( G4fs::is_directory(dirName.data()) ) {
     fIsNtupleDirectory = G4VFileManager::SetNtupleDirectoryName(dirName);
     return fIsNtupleDirectory;
  }

  G4Analysis::Warn("Directory " + dirName + " does not exists.\n"
    "Ntuples will be written in the current directory.",
    fkClass, "SetNtupleDirectoryName");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::NotifyNtupleFile(CsvNtupleDescription* ntupleDescription)
{
  // Notify not empty file
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

  return SetIsEmpty(ntupleFileName, ! ntupleDescription->GetHasFill());
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::CreateNtupleFile(
  CsvNtupleDescription* ntupleDescription)
{
  // Get ntuple file name per object (if defined)
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

  // Update file name if it is already in use
  while ( GetTFile(ntupleFileName, false) != nullptr ) {
    // the file is already in use
    auto oldName = ntupleFileName;
    auto newName = GetBaseName(oldName) + "_bis." + GetExtension(oldName);
    ntupleDescription->SetFileName(newName);

    Warn("Ntuple filename " + oldName + " is already in use.\n" +
         "It will be replaced with : " + newName,
         fkClass, "CreateNtupleFile");

    ntupleFileName = GetNtupleFileName(ntupleDescription);
  }

  // Create new ntuple file
  ntupleDescription->SetFile(CreateTFile(ntupleFileName));

  return (ntupleDescription->GetFile() != nullptr);
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::CloseNtupleFile(
  CsvNtupleDescription* ntupleDescription)
{
  // Notifying not empty file is done in G4CsvNtupleFileManager::ActionAtWrite,
  // as here we increment the cycle number and GetNtupleFileName returns a file name
  // for the next cycle version.

  // Ntuple files are registered in file manager map.
  // they will be closed with CloseFiles() calls

  ntupleDescription->GetFile().reset();

  return true;
}
