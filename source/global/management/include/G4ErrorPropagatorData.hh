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
// G4ErrorPropagatorData
//
// Class description:
//
// Utility class to provide access to mode, state, target
// and manager verbosity for the error propagation classes.

// Author: P.Arce, 2004
// --------------------------------------------------------------------

#ifndef G4ErrorPropagatorData_hh
#define G4ErrorPropagatorData_hh

#include "globals.hh"

enum G4ErrorMode
{
  G4ErrorMode_PropForwards = 1,
  G4ErrorMode_PropBackwards,
  G4ErrorMode_PropTest
};

enum G4ErrorState
{
  G4ErrorState_PreInit = 1,
  G4ErrorState_Init,
  G4ErrorState_Propagating,
  G4ErrorState_TargetCloserThanBoundary,
  G4ErrorState_StoppedAtTarget
};

enum G4ErrorStage
{
  G4ErrorStage_Inflation = 1,
  G4ErrorStage_Deflation
};

class G4ErrorTarget;

class G4ErrorPropagatorData
{
 public:
  static G4ErrorPropagatorData* GetErrorPropagatorData();
  // Singleton instance

  inline G4ErrorMode GetMode() const;
  inline void SetMode(G4ErrorMode mode);

  inline G4ErrorState GetState() const;
  inline void SetState(G4ErrorState sta);

  inline G4ErrorStage GetStage() const;
  inline void SetStage(G4ErrorStage sta);

  inline const G4ErrorTarget* GetTarget(G4bool mustExist = false) const;
  inline void SetTarget(const G4ErrorTarget* target);

  static G4int verbose();
  static void SetVerbose(G4int ver);

 private:
  G4ErrorPropagatorData() = default;
  ~G4ErrorPropagatorData();
  // constructor and destructor are private

 private:
  static G4ThreadLocal G4ErrorPropagatorData* fpInstance;

  G4ErrorMode theMode{G4ErrorMode_PropTest};

  G4ErrorState theState{G4ErrorState_PreInit};

  G4ErrorStage theStage{G4ErrorStage_Inflation};

  G4ErrorTarget* theTarget = nullptr;

  static G4ThreadLocal G4int theVerbosity;
};

#include "G4ErrorPropagatorData.icc"

#endif
