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
// $Id: G4XmlFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4XmlFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/waxml/begend"

using namespace G4Analysis;

//_____________________________________________________________________________
G4XmlFileManager::G4XmlFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fHnFile(nullptr)
{
}

//_____________________________________________________________________________
G4XmlFileManager::~G4XmlFileManager()
{}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlFileManager::OpenFile(const G4String& fileName)
{
  // Keep and locks file name
  fFileName =  fileName;
  fLockFileName = true;
  fIsOpenFile = true;

  return true;
}  
  
//_____________________________________________________________________________
G4bool G4XmlFileManager::WriteFile() 
{
  // Nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseFile()
{
  // Unlock file name
  
  fLockFileName = false;
  fIsOpenFile = false;
  return true; 
} 
   
//_____________________________________________________________________________
G4bool G4XmlFileManager::CreateHnFile()
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "histo file", GetFullFileName());
#endif
  
  // delete a previous file if it exists
  //if ( fHnFile ) delete fHnFile; 
  
  fHnFile = std::make_shared<std::ofstream>(GetFullFileName());
  if ( fHnFile->fail() ) {
    fHnFile = nullptr;
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << GetFullFileName();
    G4Exception("G4XmlFileManager::CreateHnFile()",
              "Analysis_W001", JustWarning, description);
    return false;
  }

  tools::waxml::begin(*fHnFile);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("create", "histo file", GetFullFileName());
#endif

  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseHnFile()
{
  // No file may be open if no master manager is instantiated
  // and no histograms were booked
  if ( ! fHnFile.get() ) return true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("close", "histo file", GetFullFileName());
#endif

  // close file
  tools::waxml::end(*fHnFile);
  fHnFile->close(); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("close", "histo file", GetFullFileName());
#endif


  return true; 
} 
   
//_____________________________________________________________________________
G4bool G4XmlFileManager::CreateNtupleFile(
  G4TNtupleDescription<tools::waxml::ntuple>* ntupleDescription)
{
  G4String ntupleName = ntupleDescription->fNtupleBooking.name();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("create", "ntuple file", GetNtupleFileName(ntupleName));
#endif

  auto ntupleFile = new std::ofstream(GetNtupleFileName(ntupleName));
  if ( ntupleFile->fail() ) {
    delete ntupleFile;
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " 
                << GetNtupleFileName(ntupleName);
    G4Exception("G4XmlFileManager::CreateNtupleFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }
  
  tools::waxml::begin(*ntupleFile);
  ntupleDescription->fFile = ntupleFile;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("create", "ntuple file", GetNtupleFileName(ntupleName));
#endif

  return true;
}  

//_____________________________________________________________________________
G4bool G4XmlFileManager::CloseNtupleFile(
  G4TNtupleDescription<tools::waxml::ntuple>* ntupleDescription)
{
  // Do nothing if there is no file
  if ( ! ntupleDescription->fFile ) return true;

  G4String ntupleName = ntupleDescription->fNtupleBooking.name();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()
      ->Message("close", "ntuple file", GetNtupleFileName(ntupleName));
#endif

  // close file
  tools::waxml::end(*(ntupleDescription->fFile));
  ntupleDescription->fFile->close(); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("close", "ntuple file", GetNtupleFileName(ntupleName));
#endif

  return true; 
} 
   
