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
  : fState(state),
    fFileName(""), 
    fHistoDirectoryName(""), 
    fNtupleDirectoryName(""),
    fLockFileName(false),
    fLockHistoDirectoryName(false), 
    fLockNtupleDirectoryName(false)
{
}

//_____________________________________________________________________________
G4VFileManager::~G4VFileManager()
{
}

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
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fFileName = fileName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VFileManager::SetHistoDirectoryName(const G4String& dirName) 
{
  if ( fLockHistoDirectoryName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Histo directory name as its value was already used.";
    G4Exception("G4VFileManager::SetHistoDirectoryName()",
                "Analysis_W009", JustWarning, description);
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
                "Analysis_W010", JustWarning, description);
    return false;
  }              

  fNtupleDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4String G4VFileManager::GetFullFileName() const 
{  
  G4String name(fFileName);
  // Add thread Id to a file name if MT processing
  if ( ! fState.GetIsMaster() ) {
    std::ostringstream os;
    os << G4Threading::G4GetThreadId();
    name.append("_t");
    name.append(os.str());
  }  
  // Add file extension .root if no extension is given
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }  

  return name;
}  

//_____________________________________________________________________________
G4String G4VFileManager::GetNtupleFileName(const G4String& ntupleName) const 
{  
  G4String name(fFileName);
  // Add ntupleName
  name.append("_");
  name.append(ntupleName);
  // Add thread Id to a file name if MT processing
  if ( ! fState.GetIsMaster() ) {
    std::ostringstream os;
    os << G4Threading::G4GetThreadId();
    name.append("_t");
    name.append(os.str());
  }  
  // Add file extension .xml if no extension is given
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }
  return name;
}  

//_____________________________________________________________________________
G4String G4VFileManager::GetFileType() const 
{
  G4String fileType = fState.GetType();
  fileType.toLower();
  return fileType;
}                 


