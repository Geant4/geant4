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

// Base class for the object managers.
// It handles first object ID and its lock and provides common utility methods.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4BaseAnalysisManager_h
#define G4BaseAnalysisManager_h 1

#include "G4Fcn.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include "G4AnalysisManagerState.hh"

#include <string_view>

class G4BaseAnalysisManager
{
  public:
    explicit G4BaseAnalysisManager(const G4AnalysisManagerState& state);
    G4BaseAnalysisManager() = delete;
    virtual ~G4BaseAnalysisManager() = default;

    // The ids of objects are generated automatically
    // starting from 0; with the following function it is possible to
    // change the first Id to start from other value
    G4bool SetFirstId(G4int firstId);
    void  SetLockFirstId(G4bool lockFirstId);

    // Access method
    G4int GetFirstId() const;
    G4int GetCycle() const;

  protected:
    // Methods for verbose
    G4bool IsVerbose(G4int verboseLevel) const;
    void Message(G4int level,
                 const G4String& action,
                 const G4String& objectType,
                 const G4String& objectName = "",
                 G4bool success = true) const;

    // Data members
    const G4AnalysisManagerState& fState;
    G4int    fFirstId { 0 };
    G4bool   fLockFirstId { false };

  private:
    // Static data members
    static constexpr std::string_view fkClass { "G4BaseAnalysisManager" };
};

// inline functions

inline G4bool G4BaseAnalysisManager::IsVerbose(G4int verboseLevel) const
{ return fState.IsVerbose(verboseLevel); }

inline void G4BaseAnalysisManager::Message(
  G4int level, const G4String& action, const G4String& objectType,
  const G4String& objectName, G4bool success) const
{
  fState.Message(level, action, objectType, objectName, success);
}

inline void G4BaseAnalysisManager::SetLockFirstId(G4bool lockFirstId) {
  fLockFirstId = lockFirstId;
}

inline G4int G4BaseAnalysisManager::GetFirstId() const {
  return fFirstId;
}

inline G4int G4BaseAnalysisManager::GetCycle() const {
  return fState.GetCycle();
}

#endif

