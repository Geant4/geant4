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

#include "G4VFileManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "G4Threading.hh"

//_____________________________________________________________________________
G4VFileManager::G4VFileManager(const G4AnalysisManagerState& state)
  : G4BaseFileManager(state),
    fHistoDirectoryName(""), 
    fNtupleDirectoryName(""),
    fIsOpenFile(false),
    fLockDirectoryNames(false), 
    fH1FileManager(nullptr),
    fH2FileManager(nullptr),
    fH3FileManager(nullptr),
    fP1FileManager(nullptr),
    fP2FileManager(nullptr)
{}

//_____________________________________________________________________________
G4VFileManager::~G4VFileManager()
{}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4VFileManager::SetFileName(const G4String& fileName) 
{
  // Check extension
  auto name = fileName;
  auto extension = G4Analysis::GetExtension(fileName);
  if ( extension.size() && GetFileType().size() && extension != GetFileType() ) {
    // replace extension
    name = G4Analysis::GetBaseName(fileName) + "." + GetFileType();
    G4ExceptionDescription description;
    description 
      << fileName << " file extension is not valid for " 
      << GetFileType() << " output." << G4endl
      << name << " will be used.";
    G4Exception("G4VFileManager::SetFileName()",
                "Analysis_W012", JustWarning, description);
  } 

  return G4BaseFileManager::SetFileName(name);
}

//_____________________________________________________________________________
G4bool G4VFileManager::SetHistoDirectoryName(const G4String& dirName) 
{
  if ( fLockDirectoryNames ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Histo directory name as its value was already used.";
    G4Exception("G4VFileManager::SetHistoDirectoryName()",
                "Analysis_W012", JustWarning, description);
    return false;
  }              

  fHistoDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VFileManager::SetNtupleDirectoryName(const G4String& dirName) 
{
  if ( fLockDirectoryNames ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Ntuple directory name as its value was already used.";
    G4Exception("G4VFileManager::SetNtupleDirectoryName()",
                "Analysis_W012", JustWarning, description);
    return false;
  }              

  fNtupleDirectoryName = dirName;
  return true;
}
