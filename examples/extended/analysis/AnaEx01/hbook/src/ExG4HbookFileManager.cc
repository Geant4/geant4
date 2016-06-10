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
// $Id$
//
/// \file hbook/src/ExG4HbookFileManager.cc
/// \brief Implementation of the ExG4HbookFileManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookFileManager.hh"
#include "G4AnalysisManagerState.hh"

#include <iostream>

const G4String ExG4HbookFileManager::fgkDefaultNtupleDirectoryName = "ntuple";

//_____________________________________________________________________________
ExG4HbookFileManager::ExG4HbookFileManager(const G4AnalysisManagerState& state)
 : G4VFileManager(state),
   fFile(0)
{
}

//_____________________________________________________________________________
ExG4HbookFileManager::~ExG4HbookFileManager()
{  
  delete fFile;  
}

// 
// public methods
//

//_____________________________________________________________________________
G4bool ExG4HbookFileManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;
  G4String name = GetFullFileName();

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "analysis file", name);
#endif
  
  // delete a previous file if it exists
  if ( fFile ) delete fFile; 
  
  tools::hbook::CHCDIR("//PAWC"," ");
  
  unsigned int unit = 1;
  fFile = new tools::hbook::wfile(std::cout, name, unit);
  if ( ! fFile->is_valid() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4HbookAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;       
  }

  // At this point, in HBOOK, we should have :
  //   - created a //LUN1 directory attached to the file
  //   - created a //PAWC/LUN1 in memory
  //   - be in the directory //PAWC/LUN1.

  // create an "histo" HBOOK directory both in memory and in the file :
  if ( fHistoDirectoryName != "" ) {
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    tools::hbook::CHMDIR(fHistoDirectoryName.data()," ");
    tools::hbook::CHCDIR("//LUN1"," ");
    tools::hbook::CHMDIR(fHistoDirectoryName.data()," ");
  }
  // the five upper lines could have been done with :
  //fFile->cd_home();
  //fFile->mkcd("histo");

  fLockFileName = true;
  fLockHistoDirectoryName = true;
  fLockNtupleDirectoryName = true;

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("open", "analysis file", name);
#endif
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool ExG4HbookFileManager::WriteFile() 
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("write", "file", GetFullFileName());
#endif

  // Return to //PAWC/LUN1 :
  //tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  G4bool result = fFile->write();  

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()->Message("write", "file", GetFullFileName(), result);
#endif

  return result;  
}

//_____________________________________________________________________________
G4bool ExG4HbookFileManager::CloseFile()
{
  // close file
  G4bool result = fFile->close();  
  fLockFileName = false;

  return result;
} 

//_____________________________________________________________________________
void ExG4HbookFileManager::CreateNtupleDirectory()
{
// Create an "ntuple" directory both in memory and in the file

  static G4bool isDone = false;
  
  // Do not create directory more than once
  if (isDone) return;
  
  fFile->cd_home();      //go under //PAWC/LUN1
  if ( fNtupleDirectoryName == "" )
    fFile->mkcd(fgkDefaultNtupleDirectoryName.data());
  else  
    fFile->mkcd(fNtupleDirectoryName.data());
  fLockNtupleDirectoryName = true;
  isDone = false;
}
                                     
#endif
