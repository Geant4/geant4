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
// $Id: G4VFileManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Base class for File manager.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VFileManager_h
#define G4VFileManager_h 1

#include "G4BaseFileManager.hh"
#include "globals.hh"

class G4VFileManager : public G4BaseFileManager
{
  public:
    explicit G4VFileManager(const G4AnalysisManagerState& state);
    virtual ~G4VFileManager();
   
    // Methods to manipulate files
    virtual G4bool OpenFile(const G4String& fileName) = 0;
    virtual G4bool WriteFile() = 0;
    virtual G4bool CloseFile() = 0; 
    
    // Methods for handling files and directories names
    //
    virtual G4bool SetFileName(const G4String& fileName) final;
    
    void LockHistoDirectoryName();
    void LockNtupleDirectoryName();

    G4bool SetHistoDirectoryName(const G4String& dirName);
    G4bool SetNtupleDirectoryName(const G4String& dirName); 

    G4bool IsOpenFile() const;
    G4String GetHistoDirectoryName() const;
    G4String GetNtupleDirectoryName() const;

  protected:
    // data members
    G4bool   fIsOpenFile;
    G4String fHistoDirectoryName;
    G4String fNtupleDirectoryName; 
    G4bool   fLockFileName;     
    G4bool   fLockHistoDirectoryName;     
    G4bool   fLockNtupleDirectoryName;
};

// inline functions

inline G4bool G4VFileManager::IsOpenFile() const
{ return fIsOpenFile; }

inline void G4VFileManager::LockHistoDirectoryName()
{ fLockHistoDirectoryName = true; }

inline void G4VFileManager::LockNtupleDirectoryName()
{ fLockNtupleDirectoryName = true; }


inline G4String G4VFileManager::GetHistoDirectoryName() const {
  return fHistoDirectoryName;
}  

inline G4String G4VFileManager::GetNtupleDirectoryName() const {
  return fNtupleDirectoryName;
}  
  
#endif
