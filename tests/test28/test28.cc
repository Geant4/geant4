
#include "Tst28DetectorConstruction.hh"
#include "Tst28RunAction.hh"
#include "Tst28PrimaryGeneratorAction.hh"
#include "Tst28PhysicsList.hh"
#include "Tst28SteppingAction.hh"
#include "Tst28StackingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst28DetectorConstruction);
  runManager->SetUserInitialization(new Tst28PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst28RunAction);
  runManager->SetUserAction(new Tst28PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst28StackingAction);

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}

