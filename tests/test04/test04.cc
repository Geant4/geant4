// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test04.cc,v 1.1 1999-01-08 16:34:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
//      GEANT4 - test04
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------


#include "MyDetectorConstruction.hh"
#include "MyPhysicsList.hh"
#include "MyRunAction.hh"
#include "MyPrimaryGeneratorAction.hh"
#include "MyEventAction.hh"
#include "MyStackingAction.hh"

#include "G4ios.hh"

#include "G4RunManager.hh"
#include "G4PersistencyManager.hh"

int main()
{
  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Persistency Manager
  G4PersistencyManager * persistencyManager = new G4PersistencyManager;

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  runManager->SetUserInitialization(new MyDetectorConstruction);
  runManager->SetUserInitialization(new MyPhysicsList);
  
  // UserAction classes.
  runManager->SetUserAction(new MyRunAction);
  runManager->SetUserAction(new MyPrimaryGeneratorAction);
  runManager->SetUserAction(new MyEventAction);
  runManager->SetUserAction(new MyStackingAction);

  runManager->Initialize();

  // Event loop
  G4int n_event = 10;
  runManager->BeamOn(n_event);

  delete runManager;
  return 0;
}

