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

// Base class for analysis file managers.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4BaseFileManager_h
#define G4BaseFileManager_h 1

#include "G4AnalysisManagerState.hh"
#include "globals.hh"

#include <vector>

class G4BaseFileManager
{
  public:
    explicit G4BaseFileManager(const G4AnalysisManagerState& state);
    G4BaseFileManager() = delete;
    virtual ~G4BaseFileManager() = default;

    virtual G4bool SetFileName(const G4String& fileName);
      // Set the base file name (without extension)
    virtual G4String GetFileType() const;
      // Return the manager type (starts with a lowercase letter)
    virtual G4bool HasCycles() const;
      // Return true when the output supports writing the same object in a file
      // multiple times.
      // The default implementation returns false.

    void AddFileName(const G4String& fileName);
      // For handling multiple files
      // Save file name in a vector if not yet present

    G4String GetFileName() const;
      // Return the base file name (without extension)
    G4String GetFullFileName(const G4String& baseFileName = "",
                             G4bool isPerThread = true) const;
      // Compose and return the full file name:
      // - add _tN suffix to the file base name if isPerThread
      // - add file extension if not present and available in state
    const std::vector<G4String>& GetFileNames() const;
     // Return the file names vector

    G4String GetHnFileName(const G4String& hnType,
                           const G4String& hnName) const;
      // Compose and return the histogram or profile specific file name:
      // - add _hn_hnName suffix to the file base name
      // - add file extension if not present

    G4String GetHnFileName(const G4String& fileName,
                           G4int cycle = 0) const;
      // Update Hn file name:
      // - add _vN  suffix to the base namer if cycle > 0

    G4String GetNtupleFileName(const G4String& ntupleName,
                               G4int cycle = 0) const;
      // Compose and return the ntuple specific file name:
      // - add _nt_ntupleName suffix to the file base name
      // - add _vN  suffix if cycle > 0
      // - add _tN suffix if called on thread worker
      // - add file extension if not present

    G4String GetNtupleFileName(G4int ntupleFileNumber,
                                G4int cycle = 0) const;
      // Compose and return the ntuple specific file name:
      // - add _mN suffix to the file base name
      // - add _vN  suffix if cycle > 0
      // - add file extension if not present

    G4String GetPlotFileName() const;
      // Return the file name for batch plotting output

  protected:
    // Methods for verbose
    void Message(G4int level,
                 const G4String& action,
                 const G4String& objectType,
                 const G4String& objectName = "",
                 G4bool success = true) const;

    // Data members
    const G4AnalysisManagerState& fState;
    G4String fFileName;  // to be changed in fDefaultFileName
    std::vector<G4String> fFileNames;
};

inline G4bool G4BaseFileManager::SetFileName(const G4String& fileName) {
  // CHECK if still needed in this base class
  fFileName = fileName;
  return true;
}

inline G4bool G4BaseFileManager::HasCycles() const {
  return false;
}

inline G4String G4BaseFileManager::GetFileName() const {
  return fFileName;
}

inline const std::vector<G4String>& G4BaseFileManager::GetFileNames() const {
  return fFileNames;
}

inline void G4BaseFileManager::Message(
  G4int level, const G4String& action, const G4String& objectType,
  const G4String& objectName, G4bool success) const
{
  fState.Message(level, action, objectType, objectName, success);
}

#endif
