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
// $Id: G4VFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4VFileManager.hh"
#include "G4AnalysisManagerState.hh"

#include "G4Threading.hh"

//_____________________________________________________________________________
G4VFileManager::G4VFileManager(const G4AnalysisManagerState& state)
  : G4BaseFileManager(state),
    fIsOpenFile(false),
    fHistoDirectoryName(""), 
    fNtupleDirectoryName(""),
    fLockFileName(false),
    fLockHistoDirectoryName(false), 
    fLockNtupleDirectoryName(false)
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
  if ( fLockFileName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set File name as its value was already used.";
    G4Exception("G4VFileManager::SetFileName()",
                "Analysis_W012", JustWarning, description);
    return false;
  }              

  return G4BaseFileManager::SetFileName(fileName);
}  

//_____________________________________________________________________________
G4bool G4VFileManager::SetHistoDirectoryName(const G4String& dirName) 
{
  if ( fLockHistoDirectoryName ) {
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
  if ( fLockNtupleDirectoryName ) {
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
