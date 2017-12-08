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
// $Id: G4BaseFileManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 10/09/2014  (ivana@ipno.in2p3.fr)

#include "G4BaseFileManager.hh"

#include "G4Threading.hh"

//_____________________________________________________________________________
G4BaseFileManager::G4BaseFileManager(const G4AnalysisManagerState& state)
  : fState(state),
    fFileName("")
{
}

//_____________________________________________________________________________
G4BaseFileManager::~G4BaseFileManager()
{
}

// 
// private methods
//

//_____________________________________________________________________________
G4String G4BaseFileManager::TakeOffExtension(G4String& name) const
{
  G4String extension;
  if ( name.rfind(".") != std::string::npos ) { 
    extension = name.substr(name.rfind("."));
    name = name.substr(0, name.rfind("."));
  }
  else {
    extension = ".";
    extension.append(GetFileType());
  }
  return extension;
}    

// 
// public methods
//

//_____________________________________________________________________________
G4String G4BaseFileManager::GetFullFileName(const G4String& baseFileName,
                                            G4bool isPerThread) const 
{  
  G4String name(baseFileName);
  if ( name == "" ) name = fFileName;

  // Take out file extension
  G4String extension = TakeOffExtension(name);
    
  // Add thread Id to a file name if MT processing
  if ( isPerThread && ! fState.GetIsMaster() ) {
    std::ostringstream os;
    os << G4Threading::G4GetThreadId();
    name.append("_t");
    name.append(os.str());
  }  

  // Add (back if it was present) file extension
  name.append(extension);

  return name;
}  

//_____________________________________________________________________________
G4String G4BaseFileManager::GetNtupleFileName(const G4String& ntupleName) const 
{  
  G4String name(fFileName);

  // Take out file extension
  auto extension = TakeOffExtension(name);
    
  // Add ntupleName
  name.append("_nt_");
  name.append(ntupleName);

  // Add thread Id to a file name if MT processing
  if ( ! fState.GetIsMaster() ) {
    std::ostringstream os;
    os << G4Threading::G4GetThreadId();
    name.append("_t");
    name.append(os.str());
  }  

  // Add (back if it was present) file extension
  name.append(extension);
  
  return name;
}  

//_____________________________________________________________________________
G4String G4BaseFileManager::GetNtupleFileName(G4int ntupleFileNumber) const 
{  
  G4String name(fFileName);

  // Take out file extension
  auto extension = TakeOffExtension(name);
    
  // Add _M followed by ntupleFileNumber
  std::ostringstream os;
  os << ntupleFileNumber;
  name.append("_m");
  name.append(os.str());

  // Add (back if it was present) file extension
  name.append(extension);
  
  return name;
}  

//_____________________________________________________________________________
G4String G4BaseFileManager::GetHnFileName(const G4String& hnType,
                                          const G4String& hnName) const
{
  G4String name(fFileName);

  // Take out file extension
  auto extension = TakeOffExtension(name);
 
  // Add _hnType_hnName
  name.append("_");
  name.append(hnType);
  name.append("_");
  name.append(hnName);

  // Add (back if it was present) file extension
  name.append(extension);

  return name;
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetPlotFileName() const
{
  G4String name(fFileName);

  // Take out file extension
  auto extension = TakeOffExtension(name);

  // Add .ps extension
  name.append(".ps");

  return name;
}

//_____________________________________________________________________________
G4String G4BaseFileManager::GetFileType() const 
{
  G4String fileType = fState.GetType();
  fileType.toLower();
  return fileType;
}                 


