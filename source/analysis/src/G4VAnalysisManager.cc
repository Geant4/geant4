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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4VAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4VAnalysisManager::G4VAnalysisManager(const G4String& type)
  : fVerboseLevel(0),
    fFirstHistoId(0),
    fFirstNtupleColumnId(0),
    fHistoDirectoryName(""), 
    fNtupleDirectoryName(""),
    fLockFirstHistoId(false),
    fLockFirstNtupleColumnId(false),
    fLockHistoDirectoryName(false), 
    fLockNtupleDirectoryName(false),
    fVerboseL1(type,1),
    fVerboseL2(type,2),
    fVerboseL3(type,3),
    fpVerboseL1(0),
    fpVerboseL2(0),
    fpVerboseL3(0)
{}

//_____________________________________________________________________________
G4VAnalysisManager::~G4VAnalysisManager()
{}

//_____________________________________________________________________________
void G4VAnalysisManager::SetVerboseLevel(G4int verboseLevel) 
{
  if ( verboseLevel == fVerboseLevel || verboseLevel < 0 ) return;
  
  fVerboseLevel = verboseLevel;
  
  if ( verboseLevel == 0 ) {
    fpVerboseL1 = 0;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
  }
  else if ( verboseLevel == 1 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
  }
  else if ( verboseLevel == 2 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = 0;
  }
  else {
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = &fVerboseL3;
  }
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetHistoDirectoryName(const G4String& dirName) 
{
  if ( fLockHistoDirectoryName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Histo directory name as its value was already used.";
    G4Exception("G4VAnalysisManager::SetHistoDirectoryName()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fHistoDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetNtupleDirectoryName(const G4String& dirName) 
{
  if ( fLockNtupleDirectoryName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Ntuple directory name as its value was already used.";
    G4Exception("G4VAnalysisManager::SetNtupleDirectoryName()",
                "Analysis_W010", JustWarning, description);
    return false;
  }              

  fNtupleDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstHistoId(G4int firstId) 
{
  if ( fLockFirstHistoId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstHistoId as its value was already used.";
    G4Exception("G4VAnalysisManager::SetFirstHistoId()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fFirstHistoId = firstId;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstNtupleColumnId(G4int firstId) 
{
  if ( fLockFirstNtupleColumnId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstNtupleColumnId as its value was already used.";
    G4Exception("G4VAnalysisManager::SetFirstHistoId()",
                "Analysis_W010", JustWarning, description);
    return false;
  }              

  fFirstNtupleColumnId = firstId;
  return true;
}

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetFileType() const {
  G4String fileType = fVerboseL1.GetType();
  fileType.toLower();
  return fileType;
}                 

