// $Id: exampleE03.cc,v 1.1 1999/01/07 16:05:25 gunter Exp $
#include "ExE03DetectorConstruction.hh"
#include "ExE03RunAction.hh"
#include "ExE03PrimaryGeneratorAction.hh"
#include "ExE03PhysicsList.hh"
#include "G4AssemblyCreator.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

int main(int argc,char** argv) {


  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new ExE03DetectorConstruction);
  runManager->SetUserInitialization(new ExE03PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new ExE03RunAction);
  runManager->SetUserAction(new ExE03PrimaryGeneratorAction);

  // User interactions
  // Define (G)UI for interactive mode
  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
     G4UIsession * session = new G4UIterminal;
     session->SessionStart();
     delete session;
  }
  else
  // Batch mode
  {
    G4UImanager * UI = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}







