// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisStateDependent.cc,v 1.1 1999-11-29 15:20:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4VisStateDependent.hh"

#include "G4VisManager.hh"
#include "G4StateManager.hh"

G4VisStateDependent::G4VisStateDependent (G4VisManager* pVisManager):
  fpVisManager (pVisManager) {}

G4bool G4VisStateDependent::Notify (G4ApplicationState requestedState) {
  G4StateManager* stateManager = G4StateManager::GetStateManager ();
  if(stateManager -> GetPreviousState () == EventProc &&
     requestedState == GeomClosed) {
    fpVisManager -> EndOfEvent ();
  }
  return true;
}
