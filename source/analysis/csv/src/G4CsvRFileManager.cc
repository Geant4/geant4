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
// $Id: G4CsvRFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 21/10/2014  (ivana@ipno.in2p3.fr)

#include "G4CsvRFileManager.hh"
#include "G4AnalysisManagerState.hh"

//_____________________________________________________________________________
G4CsvRFileManager::G4CsvRFileManager(const G4AnalysisManagerState& state)
 : G4BaseFileManager(state),
   fRFiles()
{
}

//_____________________________________________________________________________
G4CsvRFileManager::~G4CsvRFileManager()
{  
  for (G4int i=0; i<G4int(fRFiles.size()); ++i) { 
    delete fRFiles[i];
  }   
}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvRFileManager::OpenRFile(const G4String& fileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "read analysis file", fileName);
#endif

  // create new file
  std::ifstream* newFile = new std::ifstream(fileName);
  if ( ! newFile->is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4CsvAnalysisReader::OpenRFile()",
                "Analysis_WR001", JustWarning, description);
    return false;
  }

  // add file in a map and delete the previous file if it exists
  std::map<G4String, std::ifstream*>::iterator it
    = fRFiles.find(fileName);
  if ( it != fRFiles.end() ) { 
    delete it->second;
    it->second = newFile;
  }
  else {
    fRFiles[fileName] = newFile;
  }   

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read analysis file", fileName);
#endif

  return true;
}  
  
//_____________________________________________________________________________
std::ifstream* G4CsvRFileManager::GetRFile(const G4String& fileName) const
{ 
  std::map<G4String, std::ifstream*>::const_iterator it
    = fRFiles.find(fileName);
  if  ( it != fRFiles.end() )
    return it->second;
  else {
    return nullptr;
  }     
}
