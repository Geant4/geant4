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
// $Id: G4XmlRFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#include "G4XmlRFileManager.hh"
#include "G4AnalysisManagerState.hh"

#include "tools/waxml/begend"

//_____________________________________________________________________________
G4XmlRFileManager::G4XmlRFileManager(const G4AnalysisManagerState& state)
 : G4BaseFileManager(state),
   fReadFactory(0),
   fRFiles()
{
}

//_____________________________________________________________________________
G4XmlRFileManager::~G4XmlRFileManager()
{  
  for (G4int i=0; i<G4int(fRFiles.size()); ++i) { 
    delete fRFiles[i];
  }   
  delete fReadFactory; 
}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4XmlRFileManager::OpenRFile(const G4String& fileName)
{
  // Get full file name (add only extension)
  G4bool isPerThread = false;
  G4String name = GetFullFileName(fileName, isPerThread);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("open", "read analysis file", name);
#endif

  G4bool verbose = false;

  // create factory (if it does not yet exist)
  if ( ! fReadFactory ) {
    fReadFactory = new tools::xml::default_factory();
  }  

  // create new file
  tools::raxml* newFile 
    = new tools::raxml(*fReadFactory, G4cout, verbose);

  // clear objects 
  // (this should not be needed when starting a new raxml)
  std::vector<tools::raxml_out>& objs = newFile->objects();
  objs.clear();

  G4bool compressed = false;
  if ( ! newFile->load_file(name, compressed) ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << name;
    G4Exception("G4XmlRFileManager::OpenRFile()",
                "Analysis_WR001", JustWarning, description);
    delete newFile;
    return false;
  }

  // add file in a map and delete the previous file if it exists
  std::map<G4String, tools::raxml*>::iterator it
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
tools::raxml* G4XmlRFileManager::GetRFile(const G4String& fileName) const
{ 
  // Get full file name (add only extension)
  G4bool isPerThread = false;
  G4String name = GetFullFileName(fileName, isPerThread);

  std::map<G4String, tools::raxml*>::const_iterator it
    = fRFiles.find(name);
  if  ( it != fRFiles.end() )
    return it->second;
  else {
    return nullptr;
  }     
}
