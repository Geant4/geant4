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
// $Id: G4ErrorPropagatorManager.hh 66892 2013-01-17 10:57:59Z gunter $
//
// Class Description:
//
// This is the class manager of the GEANT4e package.
// It is the main interface for the user to define the setup and
// start the propagation. Initializes GEANT4 for the propagation.

// History:
// - Created:   Pedro Arce, February 2001
// --------------------------------------------------------------------

#ifndef G4ErrorPropagatorManager_hh
#define G4ErrorPropagatorManager_hh

#include "globals.hh"
#include "G4ErrorPropagatorData.hh"
#include "G4ErrorPropagator.hh"
#include "G4ApplicationState.hh"

class G4ErrorPropagationNavigator;
class G4ErrorRunManagerHelper;
class G4ErrorTarget;
class G4ErrorTrajState;

class G4VUserDetectorConstruction;
class G4VPhysicalVolume;
class G4VUserPhysicsList;
class G4UserRunAction;
class G4UserEventAction;
class G4UserStackingAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4Mag_UsualEqRhs;
class G4Track;

class G4ErrorPropagatorManager
{

 public:  // with description

  G4ErrorPropagatorManager();
    // Initialise data to 0. Starts the G4ErrorRunManagerHelper and
    // G4ErrorPropagationNavigator.

  ~G4ErrorPropagatorManager();

  static G4ErrorPropagatorManager* GetErrorPropagatorManager();
    // Get only instance of G4ErrorPropagatorManager. If it does not exists,
    // creates it

  void EventTermination();
    // Set state to G4ErrorState_Init

  void RunTermination();
    // Set state to G4ErrorState_Init and invoke
    // G4ErrorRunManagerHelper::RunTermination()

  void InitGeant4e();
    // Initializes Geant4 and Geant4e

  void InitTrackPropagation();
    // Set the propagator step number to 0 and the G4ErrorState to Propagating
 
  G4bool InitFieldForBackwards();
    // Creates the G4ErrorMag_UsualEqRhs, that will control backwards tracking

  G4int Propagate( G4ErrorTrajState* currentTS, const G4ErrorTarget* target, G4ErrorMode mode = G4ErrorMode_PropForwards );
    // Inits track propagation, invokes G4ErrorPropagator::Propagate and
    // terminates "event"

  G4int PropagateOneStep( G4ErrorTrajState* currentTS,
                          G4ErrorMode mode = G4ErrorMode_PropForwards );
    // Invokes G4ErrorPropagator::PropagateOneStep

  G4bool CloseGeometry();
    // Close Geant4 geometry

  void SetUserInitialization(G4VUserDetectorConstruction* userInit);
    // Invokes G4ErrorRunManagerHelper to construct detector and set
    // world volume
  void SetUserInitialization(G4VPhysicalVolume* userInit);
    // Invokes G4ErrorRunManagerHelper to  set world volume
  void SetUserInitialization(G4VUserPhysicsList* userInit);
    // Invokes G4ErrorRunManagerHelper to initialize physics

  void SetUserAction(G4UserTrackingAction* userAction);
    // Invokes G4EventManager to set a G4UserTrackingAction
  void SetUserAction(G4UserSteppingAction* userAction);
    // Invokes G4EventManager to set a G4UserSteppingAction

  G4String PrintG4ErrorState();
  G4String PrintG4ErrorState( G4ErrorState state );
    // Print Geant4e state

  G4String PrintG4State();
  G4String PrintG4State( G4ApplicationState state );
    // Print Geant4 state
 
  // Set and Get methods 

  G4ErrorRunManagerHelper* GetErrorRunManagerHelper() const
    { return theG4ErrorRunManagerHelper; }

  void SetSteppingManagerVerboseLevel();

  G4ErrorPropagationNavigator* GetErrorPropagationNavigator() const
    { return theG4ErrorPropagationNavigator; }

  G4ErrorPropagator* GetPropagator() const
    { return thePropagator; }

 private:

  void StartG4ErrorRunManagerHelper();
    // Create a G4ErrorRunManagerHelper if it does not exist and set to it
    // the G4ErrorPhysicsList
 
  void StartNavigator();
    // create a G4ErrorPropagationNavigator
 
 private:

  static G4ThreadLocal G4ErrorPropagatorManager* theG4ErrorPropagatorManager;

  G4ErrorRunManagerHelper* theG4ErrorRunManagerHelper;

  G4ErrorPropagator* thePropagator;

  G4Mag_UsualEqRhs* theEquationOfMotion;

  G4ErrorPropagationNavigator* theG4ErrorPropagationNavigator;

};

#endif
