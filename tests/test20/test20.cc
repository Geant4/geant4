// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test20.cc,v 1.2 2001-05-25 10:20:27 pia Exp $
//

#include "Tst20DetectorConstruction.hh"
#include "Tst20RunAction.hh"
#include "Tst20PrimaryGeneratorAction.hh"
#include "Tst20PhysicsList.hh"
#include "Tst20SteppingAction.hh"
#include "Tst20TrackingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  HepRandom::setTheEngine(new RanecuEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  Tst20DetectorConstruction* detector;
  detector = new Tst20DetectorConstruction;

  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst20PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst20PrimaryGeneratorAction(detector));
  Tst20RunAction* runaction = new Tst20RunAction;
  runManager->SetUserAction(runaction);

  runManager->SetUserAction(new Tst20SteppingAction);
  //runManager->SetUserAction(new Tst20TrackingAction);

  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    UImanager->ApplyCommand("/control/execute prerunTst20.mac");
    session->SessionStart();
    delete session;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}
