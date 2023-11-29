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

// The state of the analysis manager instance.
//
// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisManagerState_h
#define G4AnalysisManagerState_h 1

#include "G4AnalysisVerbose.hh"
#include "G4Threading.hh"
#include "globals.hh"

#include <string_view>

class G4AnalysisManagerState
{
  // Only G4VAnalysisManager can change the state
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;
  friend class G4ParameterManager;

  public:
    G4AnalysisManagerState(G4String type, G4bool isMaster);
    // disabled constructors, operators
    G4AnalysisManagerState() = delete;
    G4AnalysisManagerState(const G4AnalysisManagerState&) = delete;
    G4AnalysisManagerState& operator=(const G4AnalysisManagerState&) = delete;

    // Methods
    void Message([[maybe_unused]] G4int level,
                 [[maybe_unused]] const G4String& action,
                 [[maybe_unused]] const G4String& objectType,
                 [[maybe_unused]] const G4String& objectName = "",
                 [[maybe_unused]] G4bool success = true) const;
    void IncrementCycle();
    void ResetCycle();

    // get methods
    G4String GetType() const;
    G4bool   GetIsMaster() const;
    G4int    GetThreadId() const;
    G4bool   GetIsActivation() const;
    G4int    GetVerboseLevel() const;
    G4bool   IsVerbose(G4int verboseLevel) const;
    G4int    GetCompressionLevel() const;
    G4int    GetCycle() const;

  private:
    // set methods
    // (hidden from all clients except for G4VAnalysisManager friend)
    void SetIsActivation(G4bool isActivation);
    void SetVerboseLevel(G4int verboseLevel);
    void SetCompressionLevel(G4int level);

    // Static data members
    static constexpr std::string_view fkClass { "G4AnalysisManagerState" };

    // Data members
    G4String fType;
    G4bool   fIsMaster;
    G4int    fThreadId { G4Threading::SEQUENTIAL_ID };
    G4bool   fIsActivation { false };
    G4int    fVerboseLevel { 0 };
    G4int    fCompressionLevel { 1 };
    G4int    fCycle { 0 };
    G4AnalysisVerbose fVerbose;
};

// inline functions

inline void G4AnalysisManagerState::SetIsActivation(G4bool isActivation)
{ fIsActivation = isActivation; }

inline void G4AnalysisManagerState::SetCompressionLevel(G4int level)
{ fCompressionLevel = level; }

inline void G4AnalysisManagerState::IncrementCycle()
{ ++fCycle; }

inline void G4AnalysisManagerState::ResetCycle()
{ fCycle = 0; }

inline G4String G4AnalysisManagerState::GetType() const
{ return fType; }

inline G4bool  G4AnalysisManagerState::GetIsMaster() const
{ return fIsMaster; }

inline G4int  G4AnalysisManagerState::GetThreadId() const
{ return fThreadId; }

inline G4bool  G4AnalysisManagerState::GetIsActivation() const
{ return fIsActivation; }

inline G4int   G4AnalysisManagerState::GetVerboseLevel() const
{ return fVerboseLevel; }

inline G4bool  G4AnalysisManagerState::IsVerbose(G4int verboseLevel) const
{ return fVerboseLevel == verboseLevel; }

inline G4int  G4AnalysisManagerState::GetCompressionLevel() const
{ return fCompressionLevel; }

inline G4int G4AnalysisManagerState::GetCycle() const
{ return fCycle; }

#endif
