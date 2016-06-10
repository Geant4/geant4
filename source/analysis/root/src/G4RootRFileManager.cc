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
// $Id: G4RootRFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#include "G4RootRFileManager.hh"
#include "G4AnalysisManagerState.hh"

#include "tools/rroot/file"
#include <tools/zlib>

#include <iostream>
#include <cstdio>

//_____________________________________________________________________________
G4RootRFileManager::G4RootRFileManager(const G4AnalysisManagerState& state)
 : G4BaseFileManager(state),
   fRFiles()
{
}

//_____________________________________________________________________________
G4RootRFileManager::~G4RootRFileManager()
{  
  for (G4int i=0; i<G4int(fRFiles.size()); ++i) { 
    delete fRFiles[i];
  }    
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4RootRFileManager::OpenRFile(const G4String& fileName,
                                     G4bool isPerThread)
{
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "read analysis file", name);
#endif

  // create new file
  tools::rroot::file* newFile = new tools::rroot::file(G4cout, name);
  newFile->add_unziper('Z',tools::decompress_buffer);
  
  if ( ! newFile->is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << name;
    G4Exception("G4RootRFileManager::OpenFile()",
                "Analysis_WR001", JustWarning, description);
    delete newFile;
    return false;
  }

  // add file in a map and delete the previous file if it exists
  std::map<G4String, tools::rroot::file*>::iterator it
    = fRFiles.find(name);
  if ( it != fRFiles.end() ) { 
    delete it->second;
    it->second = newFile;
  }
  else {
    fRFiles[name] = newFile;
  }    

#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read analysis file", name);
#endif

  return true;
}  
  
//_____________________________________________________________________________
tools::rroot::file* G4RootRFileManager::GetRFile(const G4String& fileName,
                                                 G4bool isPerThread) const
{ 
  // Get full file name
  G4String name = GetFullFileName(fileName, isPerThread);

  std::map<G4String, tools::rroot::file*>::const_iterator it
    = fRFiles.find(name);
  if  ( it != fRFiles.end() )
    return it->second;
  else {
    return nullptr;
  }     
}
