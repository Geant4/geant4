// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StateManager.hh,v 1.1 1999-01-07 16:09:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4StateManager ----------------
//             by Gabriele Cosmo, November 1996
// Class responsible for handling and updating the running state
// of the Geant4 application during its different phases.
// The class is a singleton, it can be accessed via the public
// method G4StateManager::GetStateManager().
// ------------------------------------------------------------

#ifndef G4StateManager_h
#define G4StateManager_h 1

#include <rw/tpordvec.h>
#include "globals.hh"
#include "G4ApplicationState.hh"
#include "G4VStateDependent.hh"


class G4StateManager {

public:
  static G4StateManager* GetStateManager();

protected:
  G4StateManager();

public:
  ~G4StateManager();

public:
  G4ApplicationState GetCurrentState();
  G4ApplicationState GetPreviousState();
  G4bool SetNewState(G4ApplicationState requestedState);
  G4bool RegisterDependent(G4VStateDependent* aDependent);
  G4bool DeregisterDependent(G4VStateDependent* aDependent);
  G4VStateDependent* RemoveDependent(const G4VStateDependent* aDependent);
  G4String GetStateString(G4ApplicationState aState);
  void Pause();
  void Pause(char* msg);
  void Pause(G4String msg);

private:
  G4StateManager(const G4StateManager &right);
  G4StateManager& operator=(const G4StateManager &right);
  G4int operator==(const G4StateManager &right) const;
  G4int operator!=(const G4StateManager &right) const;

private:
  static G4StateManager* theStateManager;
  G4ApplicationState theCurrentState;
  G4ApplicationState thePreviousState;
  RWTPtrOrderedVector<G4VStateDependent> theDependentsList;

};

#endif
