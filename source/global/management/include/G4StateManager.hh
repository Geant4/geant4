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
// $Id: G4StateManager.hh 108486 2018-02-15 14:47:25Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//
//       ---------------- G4StateManager ----------------
//
// Authors: G.Cosmo, M.Asai - November 1996
//
// -------------------------------------------------------------
//
// Class Description:
//
// Class responsible for handling and updating the running state
// of the Geant4 application during its different phases.
// The class is a singleton, it can be accessed via the public
// method G4StateManager::GetStateManager().
//
// States of Geant4 are defined in G4ApplicationState.

// -------------------------------------------------------------

#ifndef G4StateManager_h
#define G4StateManager_h 1

#include <vector>
#include "globals.hh"
#include "G4ApplicationState.hh"
#include "G4VStateDependent.hh"
#include "G4VExceptionHandler.hh"

class G4StateManager
{

public: // with description

  static G4StateManager* GetStateManager();
    // The G4StateManager class is a singleton class and the pointer
    // to the only one existing object can be obtained by this static 
    // method.

protected:

  G4StateManager();

public:

  ~G4StateManager();

public: // with description

  G4ApplicationState GetCurrentState() const;
    // Returns the current state
  G4ApplicationState GetPreviousState() const;
    // Returns the previous state
  G4bool SetNewState(G4ApplicationState requestedState);
    // Set Geant4 to a new state.
    // In case the request is irregal, false will be returned
    // and the state of Geant4 will not be changed.
  G4bool SetNewState(G4ApplicationState requestedState, const char* msg);
    // Set Geant4 to a new state.
    // In case the request is irregal, false will be returned
    // and the state of Geant4 will not be changed.
    // "msg" is information associating to this state change
  G4bool RegisterDependent(G4VStateDependent* aDependent,G4bool bottom=false);
    // Register a concrete class of G4VStateDependent.
    // Registered concrete classes will be notified via
    // G4VStateDependent::Notify() method when the state of Geant4 changes.
    // False will be returned if registration fails.
  G4bool DeregisterDependent(G4VStateDependent* aDependent);
    // Remove the registration.
    // False will be returned if aDependent has not been registered.
  G4VStateDependent* RemoveDependent(const G4VStateDependent* aDependent);
    // Remove the registration.
    // Removed pointer is returned.
  G4String GetStateString(G4ApplicationState aState) const;
    // Utility method which returns a string of the state name.

public:

  inline void SetSuppressAbortion(G4int i);
  inline G4int GetSuppressAbortion() const;
  inline const char* GetMessage() const;
  inline void SetExceptionHandler(G4VExceptionHandler* eh);
  inline G4VExceptionHandler* GetExceptionHandler() const;
  static void SetVerboseLevel(G4int val);

  //void Pause();
  //void Pause(const char* msg);
  //void Pause(G4String msg);
  ////  G4UIsession::pauseSession() will be invoked. The argument string "msg"
  //// will be used as a prompt characters if the session is non-graphical.
  ////  This method can be invoked by any user action class during the event
  //// loop. After the user's interactions, control goes back to the caller.

private:

  G4StateManager(const G4StateManager &right);
  G4StateManager& operator=(const G4StateManager &right);
  G4int operator==(const G4StateManager &right) const;
  G4int operator!=(const G4StateManager &right) const;

private:

  static G4ThreadLocal G4StateManager* theStateManager;
  G4ApplicationState theCurrentState;
  G4ApplicationState thePreviousState;
  std::vector<G4VStateDependent*> theDependentsList;
  G4VStateDependent* theBottomDependent;
  G4int suppressAbortion;
  const char* msgptr;
  G4VExceptionHandler* exceptionHandler;
  static G4int verboseLevel;
};

#include "G4StateManager.icc"

#endif
