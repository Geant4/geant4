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

//_____________________________________________________________________________
G4CsvFileManager::~G4CsvFileManager()
{}

// 
// private methods
//

//_____________________________________________________________________________
G4String G4CsvFileManager::GetNtupleFileName(CsvNtupleDescription* ntupleDescription)
{
  // get ntuple file name
  auto ntupleFileName = ntupleDescription->fFileName;
  if ( ntupleFileName.size() ) {
    // update filename per object per thread
    ntupleFileName = GetTnFileName(ntupleFileName, GetFileType());
  } else {
    // compose ntuple file name from the default file name
    ntupleFileName = GetNtupleFileName(ntupleDescription->fNtupleBooking.name());
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
    file = nullptr;
    G4ExceptionDescription description;
    description << "      " << "Cannot create file " << fileName;
    G4Exception("G4CsvFileManager::CreateFileImpl()",
              "Analysis_W001", JustWarning, description);
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
G4bool G4CsvFileManager::CreateNtupleFile(
  CsvNtupleDescription* ntupleDescription)
{
  // set description file name so that we can properly save to directories
  auto path = GetNtupleDirectoryName();
  if (!path.empty()) {
      path.append("/");
    }
  ntupleDescription->fFileName = path+fFileName+"_nt_"
                                  +ntupleDescription->fNtupleBooking.name();

  // get ntuple file name per object (if defined)
  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4()->Message("create", "ntuple file", ntupleFileName);
  }
#endif

  // update file name if it is already in use
  while ( GetTFile(ntupleFileName, false) ) {
    // the file is already in use
    auto oldName = ntupleDescription->fFileName;
    auto newName = GetBaseName(oldName) + "_bis." + GetExtension(oldName);
    ntupleDescription->fFileName = newName;

    G4ExceptionDescription description;
    description
      << "Ntuple filename " << oldName << " is already in use." << G4endl
      << "It will be replaced with : " << newName;
    G4Exception("G4CsvFileManager::CreateFileImpl()",
                "Analysis_W001", JustWarning, description);

    ntupleFileName = GetNtupleFileName(ntupleDescription);
  }

  ntupleDescription->fFile = CreateTFile(ntupleFileName);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    fState.GetVerboseL2()->Message("create", "ntuple file", ntupleFileName);
  }
#endif

  return (ntupleDescription->fFile != nullptr);
}

//_____________________________________________________________________________
G4bool G4CsvFileManager::CloseNtupleFile(
  CsvNtupleDescription* ntupleDescription)
{
  // Do nothing if there is no file
  if ( ! ntupleDescription->fFile ) return true;

  auto finalResult = true;

  auto ntupleFileName = GetNtupleFileName(ntupleDescription);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    fState.GetVerboseL4() ->Message("close", "ntuple file", ntupleFileName);
  }
#endif

  // close file
  auto result = CloseTFile(ntupleFileName);
  finalResult = result && finalResult;

  // Notify not empty file
  result = SetIsEmpty(ntupleFileName, ! ntupleDescription->fHasFill);
  finalResult = result && finalResult;

  // Reset file info in ntuple description
  ntupleDescription->fFile.reset();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) {
    fState.GetVerboseL2()->Message("close", "ntuple file", ntupleFileName);
  }
#endif

  return finalResult;
}
