// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test15.cc,v 1.2 1999-12-15 14:54:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst15DetectorConstruction.hh"
#include "Tst15RunAction.hh"
#include "Tst15PrimaryGeneratorAction.hh"
#include "Tst15PhysicsList.hh"
#include "Tst15SteppingAction.hh"

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
  runManager->SetUserInitialization(new Tst15DetectorConstruction);
  runManager->SetUserInitialization(new Tst15PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst15RunAction);
  runManager->SetUserAction(new Tst15PrimaryGeneratorAction);
  //runManager->SetUserAction(new Tst15SteppingAction);

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

