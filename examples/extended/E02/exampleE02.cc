// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleE02.cc,v 1.1 1999/01/07 16:05:21 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//
//
// --------------------------------------------------------------
//      GEANT4 - exampleE02
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------


#include "ExE02DetectorConstruction.hh"
#include "ExE02PhysicsList.hh"
#include "ExE02RunAction.hh"
#include "ExE02PrimaryGeneratorAction.hh"
#include "ExE02EventAction.hh"
#include "ExE02StackingAction.hh"

#include "G4ios.hh"

#include "G4RunManager.hh"
#include "G4PersistencyManager.hh"

int main()
{
  // Persistency Manager
  G4PersistencyManager * persistencyManager =
                          G4PersistencyManager::GetPersistencyManager();

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  runManager->SetUserInitialization(new ExE02DetectorConstruction);
  runManager->SetUserInitialization(new ExE02PhysicsList);
  
  // UserAction classes.
  runManager->SetUserAction(new ExE02RunAction);
  runManager->SetUserAction(new ExE02PrimaryGeneratorAction);
  runManager->SetUserAction(new ExE02EventAction);
  runManager->SetUserAction(new ExE02StackingAction);

  runManager->Initialize();

  // Event loop
  G4int n_event = 10;
  runManager->BeamOn(n_event);

  delete runManager;
  return 0;
}

