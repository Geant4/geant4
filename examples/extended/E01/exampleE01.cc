// $Id: exampleE01.cc,v 1.1 1999/01/07 16:05:17 gunter Exp $

#include "ExE01DetectorConstruction.hh"
#include "ExE01RunAction.hh"
#include "ExE01PrimaryGeneratorAction.hh"
#include "ExE01PhysicsList.hh"
#include "ExE01SteppingAction.hh"

#include "G4RunManager.hh"

int main(int ,char** ) {
  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new ExE01DetectorConstruction);
  runManager->SetUserInitialization(new ExE01PhysicsList); 

  // UserAction classes.
  runManager->SetUserAction(new ExE01RunAction);
  runManager->SetUserAction(new ExE01PrimaryGeneratorAction);
  runManager->SetUserAction(new ExE01SteppingAction);
 
  // Event loop
  G4int n_event = 1000;
  runManager->Initialize();
  runManager->BeamOn(n_event);

  delete runManager;
  return 0;
}
