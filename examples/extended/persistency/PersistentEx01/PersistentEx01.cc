// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersistentEx01.cc,v 1.4.2.1 1999/12/07 20:47:20 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
//
// --------------------------------------------------------------
//      GEANT4 - PersistentEx01
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------


#include "PersEx01DetectorConstruction.hh"
#include "PersEx01PhysicsList.hh"
#include "PersEx01RunAction.hh"
#include "PersEx01PrimaryGeneratorAction.hh"
#include "PersEx01EventAction.hh"
#include "PersEx01StackingAction.hh"

#include "G4ios.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4PersistencyManager.hh"

int main(int argc, char** argv)
{
  // Persistency Manager
  G4PersistencyManager * persistencyManager =
                          G4PersistencyManager::GetPersistencyManager();

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  runManager->SetUserInitialization(new PersEx01DetectorConstruction);
  runManager->SetUserInitialization(new PersEx01PhysicsList);
  
  // UserAction classes.
  runManager->SetUserAction(new PersEx01RunAction);
  runManager->SetUserAction(new PersEx01PrimaryGeneratorAction);
  runManager->SetUserAction(new PersEx01EventAction);
  runManager->SetUserAction(new PersEx01StackingAction);

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if(argc==1)
  // Define (G)UI terminal for interactive mode
  {
    G4UIsession* session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerunPEx01.mac");
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  // job termination
  delete runManager;
  delete persistencyManager;

  return 0;
}

