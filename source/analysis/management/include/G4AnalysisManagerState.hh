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
// $Id: G4AnalysisManagerState.hh 66310 2012-12-17 11:56:35Z ihrivnac $

// The state of the analysis manager instance.
//
// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#ifndef G4AnalysisManagerState_h
#define G4AnalysisManagerState_h 1

#include "globals.hh"
#include "G4AnalysisVerbose.hh" 
#include "G4Threading.hh"

class G4AnalysisManagerState
{
  // Only G4VAnalysisManager can change the state
  friend class G4VAnalysisManager;
  friend class G4VAnalysisReader;
  friend class G4ParameterManager;

  public: 
    G4AnalysisManagerState(const G4String& type, G4bool isMaster);

    // get methods
    G4String GetType() const;
    G4bool   GetIsMaster() const;
    G4bool   GetIsActivation() const;
    G4int    GetVerboseLevel() const;
    const G4AnalysisVerbose* GetVerboseL1() const;
    const G4AnalysisVerbose* GetVerboseL2() const;
    const G4AnalysisVerbose* GetVerboseL3() const;
    const G4AnalysisVerbose* GetVerboseL4() const;
    G4int    GetCompressionLevel() const;

  private:
    // disabled constructors, operators  
    G4AnalysisManagerState(); 
    G4AnalysisManagerState(const G4AnalysisManagerState&); 
    G4AnalysisManagerState& operator=(const G4AnalysisManagerState&); 

    // set methods
    // (hidden from all clients except for G4VAnalysisManager friend)
    void SetIsActivation(G4bool isActivation);
    void SetVerboseLevel(G4int verboseLevel);
    void SetCompressionLevel(G4int level);

    // data members
    G4String fType;
    G4bool   fIsMaster;
    G4bool   fIsActivation;
    G4int    fVerboseLevel;
    G4int    fCompressionLevel;
    G4AnalysisVerbose  fVerboseL1;
    G4AnalysisVerbose  fVerboseL2;
    G4AnalysisVerbose  fVerboseL3;
    G4AnalysisVerbose  fVerboseL4;
    G4AnalysisVerbose* fpVerboseL1;
    G4AnalysisVerbose* fpVerboseL2;
    G4AnalysisVerbose* fpVerboseL3;
    G4AnalysisVerbose* fpVerboseL4;
};

// inline functions

inline void G4AnalysisManagerState::SetIsActivation(G4bool isActivation)
{ fIsActivation = isActivation; }

inline void G4AnalysisManagerState::SetCompressionLevel(G4int level)
{ fCompressionLevel = level; }

inline G4String G4AnalysisManagerState::GetType() const
{ return fType; }

inline G4bool  G4AnalysisManagerState::GetIsMaster() const
{ return fIsMaster; }

inline G4bool  G4AnalysisManagerState::GetIsActivation() const
{ return fIsActivation; }

inline G4int   G4AnalysisManagerState::GetVerboseLevel() const
{ return fVerboseLevel; }

inline const G4AnalysisVerbose* G4AnalysisManagerState::GetVerboseL1() const
{ return fpVerboseL1; }

inline const G4AnalysisVerbose* G4AnalysisManagerState::GetVerboseL2() const
{ return fpVerboseL2; }

inline const G4AnalysisVerbose* G4AnalysisManagerState::GetVerboseL3() const
{ return fpVerboseL3; }

inline const G4AnalysisVerbose* G4AnalysisManagerState::GetVerboseL4() const
{ return fpVerboseL4; }

inline G4int  G4AnalysisManagerState::GetCompressionLevel() const
{ return fCompressionLevel; }

#endif  

