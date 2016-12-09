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
// $Id: G4BaseFileManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Base class for analysis file managers.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4BaseFileManager_h
#define G4BaseFileManager_h 1

#include "G4AnalysisManagerState.hh"
#include "globals.hh"

class G4BaseFileManager
{
  public:
    explicit G4BaseFileManager(const G4AnalysisManagerState& state);
    virtual ~G4BaseFileManager();

    virtual G4bool SetFileName(const G4String& fileName);
      // Set the base file name (without extension)
      // If the current file name is already in use
      // setting is not performed and false is returned
    
    G4String GetFileName() const;
      // Return the base file name (without extension)

    G4String GetFullFileName(const G4String& baseFileName = "",
                             G4bool isPerThread = true) const;
      // Compose and return the full file name:
      // - add _tN suffix to the file base name if isPerThread
      // - add file extension if not present

    G4String GetHnFileName(const G4String& hnType, 
                           const G4String& hnName) const;
      // Compose and return the histogram or profile specific file name:
      // - add _hn_hnName suffix to the file base name
      // - add file extension if not present

    G4String GetNtupleFileName(const G4String& ntupleName) const;
      // Compose and return the ntuple specific file name:
      // - add _nt_ntupleName suffix to the file base name
      // - add _tN suffix if called on thread worker
      // - add file extension if not present
    
    G4String GetNtupleFileName(G4int ntupleFileNumber) const;
      // Compose and return the ntuple specific file name:
      // - add _mN suffix to the file base name 
      // - add file extension if not present
    
    G4String GetPlotFileName() const;
      // Return the file name for batch plotting output

    G4String GetFileType() const;                 
     // Return the manager file type (starts with a lowercase letter)

  protected:
    // utility function
    G4String TakeOffExtension(G4String& name) const;
  
    // data members
    const G4AnalysisManagerState& fState;
    G4String fFileName;
};

inline G4bool G4BaseFileManager::SetFileName(const G4String& fileName) {
  fFileName = fileName;
  return true;
}  

inline G4String G4BaseFileManager::GetFileName() const {
  return fFileName;
}  
  
#endif

