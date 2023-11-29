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

#include "G4XmlFileManager.hh"
#include "G4XmlHnFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/waxml/begend"

using namespace G4Analysis;
using namespace tools;

//_____________________________________________________________________________
G4XmlFileManager::G4XmlFileManager(const G4AnalysisManagerState& state)
 : G4VTFileManager<std::ofstream>(state)
{
  // Create helpers defined in the base class
  fH1FileManager = std::make_shared<G4XmlHnFileManager<histo::h1d>>(this);
  fH2FileManager = std::make_shared<G4XmlHnFileManager<histo::h2d>>(this);
  fH3FileManager = std::make_shared<G4XmlHnFileManager<histo::h3d>>(this);
  fP1FileManager = std::make_shared<G4XmlHnFileManager<histo::p1d>>(this);
  fP2FileManager = std::make_shared<G4XmlHnFileManager<histo::p2d>>(this);
}

//
// private methods
//

//_____________________________________________________________________________
G4String G4XmlFileManager::GetNtupleFileName(XmlNtupleDescription* ntupleDescription)
{
  // get ntuple file name
  auto ntupleFileName = ntupleDescription->GetFileName();
  if (ntupleFileName.size() != 0u) {
    // update filename per object per thread
    ntupleFileName =
      GetTnFileName(ntupleFileName, GetFileType()/*, ntupleDescription->GetCycle()*/);
  }
  else {
    // compose ntuple file name from the default file name
    ntupleFileName = GetNtupleFileName(
      ntupleDescription->GetNtupleBooking().name()/*, ntupleDescription->GetCycle()*/);
  }

  return ntupleFileName;
}

//
// protected methods
//

//_____________________________________________________________________________
std::shared_ptr<std::ofstream> G4XmlFileManager::CreateFileImpl(const G4String& fileName)
{
  std::shared_ptr<std::ofstream> file = std::make_shared<std::ofstream>(fileName);
  if ( file->fail() ) {
    Warn(G4String("Cannot create file ") + fileName, fkClass, "CreateFileImpl");
    return nullptr;
  }

  waxml::begin(*file);
  return file;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::WriteFileImpl(std::shared_ptr<std::ofstream> /*file*/)
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseFileImpl(std::shared_ptr<std::ofstream> file)
{
  if ( ! file ) return false;

  // close file
  waxml::end(*file);
  file->close();

  return true;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlFileManager::OpenFile(const G4String& fileName)
{
  // Keep and locks file name
  fFileName =  fileName;
  auto name = GetFullFileName(fFileName);

  if ( fFile ) {
    Warn(G4String("File ") + fileName + " already exists.", fkClass, "OpenFile");
    fFile.reset();
  }

  // Create histograms file (on master)
  if ( fState.GetIsMaster() ) {
    // Create file (and save in in the file map (on master only)
    fFile = CreateTFile(name);
    if ( ! fFile) {
      Warn(G4String("Failed to create file") + fileName, fkClass, "OpenFile");
      return false;
    }
  }

  fIsOpenFile = true;

  return true;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CreateNtupleFile(
  XmlNtupleDescription* ntupleDescription)
{
  // Get ntuple file name per object (if defined)
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

  // Create ntuple file name if it does not yet exist
  auto ntupleFile = GetTFile(ntupleFileName, false);
  if ( ntupleFile == nullptr) {
    ntupleFile = CreateTFile(ntupleFileName);
  }

  ntupleDescription->SetFile(ntupleFile);

  return (ntupleDescription->GetFile() != nullptr);
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseNtupleFile(
  XmlNtupleDescription* ntupleDescription)
{
  // Notify not empty file
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);
  auto result = SetIsEmpty(ntupleFileName, ! ntupleDescription->GetHasFill());

  // Ntuple files are registered in file manager map.
  // they will be closed with CloseFiles() calls
  ntupleDescription->GetFile().reset();

  return result;
}
