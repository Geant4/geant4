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
//
// $Id: G4ErrorPropagator.hh 68050 2013-03-13 14:36:39Z gcosmo $
//
//
// Class Description:
//
//  Manages the propagation of tracks. Creates a G4Track, asks to
//  propagate it and takes also care to propagate the errors.
//  Stops the track when GEANT4 stops it or a G4ErrorTarget is reached.

// History:
// - Created:   P. Arce
// --------------------------------------------------------------------

#ifndef G4ErrorPropagator_hh
#define G4ErrorPropagator_hh

#include "globals.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4SteppingManager.hh"

class G4eErrorMatrix;
class G4Track;
class G4ErrorTrajState;
class G4ErrorFreeTrajState;
class G4ErrorTarget;
#include "globals.hh"

class G4ErrorPropagator 
{
 public:  // with description

  G4ErrorPropagator();
  ~G4ErrorPropagator(){}

  G4Track* InitG4Track( G4ErrorTrajState& initialTS );
    // Creates a G4Track from a G4ErrorTrajState

  G4int Propagate( G4ErrorTrajState* currentTS,
                   const G4ErrorTarget* target,
                   G4ErrorMode mode = G4ErrorMode_PropForwards);
    // Steers the GEANT4 propagation of a track:
    // the particle will be extrapolated until the Target is reached.
    // The final G4Track parameters will be passed to theFinalTrajState

  G4int PropagateOneStep( G4ErrorTrajState* currentTS );
    // Propagates a G4Track by one step, and then returns control to the user

  G4int MakeOneStep( G4ErrorFreeTrajState* currentTS_FREE );
    // Advance one step

  G4ErrorFreeTrajState* InitFreeTrajState( G4ErrorTrajState* currentTS );
    // Creates theCurrentTS_FREE (transforms the user G4ErrorSurfaceTrajState
    // or copies the G4ErrorFreeTrajState)

  void GetFinalTrajState( G4ErrorTrajState* currentTS, G4ErrorFreeTrajState* currentTS_FREE, const G4ErrorTarget* target );
    // After steps are done, convert the G4ErrorFreeTrajState used for error
    // propagation to the class of origin (G4ErrorFreeTrajState or
    // G4eTrajStatOnSurface)

  void InvokePreUserTrackingAction( G4Track* fpTrack );
    // Invoke the G4UserTrackingAction::PreUserTrackingAction
  void InvokePostUserTrackingAction( G4Track* fpTrack );
    // Invoke the G4UserTrackingAction::PostUserTrackingAction

  G4bool CheckIfLastStep( G4Track* aTrack );
    // Check if it is the last step for error propagation:
    //  - G4ErrorState is G4ErrorState_StoppedAtTarget
    //  - Track is OutOfWorld
    //  - G4TrackStatus is fStopAndKill

  // Get and Set methods 

  const G4ErrorTrajState* GetInitialTrajState() const
    { return theInitialTrajState; }

  G4double GetStepLength() const
    { return theStepLength; }

  void SetStepLength( const G4double sl )
    { theStepLength = sl; }

  void SetStepN( const G4int sn )
    { theStepN = sn; }

 private:

  G4int MakeSteps( G4ErrorFreeTrajState* currentTS_FREE );
    // Advance steps until target is reached

 private:

  G4double theStepLength;

  G4ErrorTrajState* theInitialTrajState;

  G4int theStepN;

  G4Track* theG4Track;

  G4SteppingManager* fpSteppingManager;

  G4int verbose;

  G4bool thePropIsInitialized;

};

#endif
