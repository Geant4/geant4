// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test16.cc,v 1.1 1999-11-18 14:58:15 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst16DetectorConstruction.hh"
#include "Tst16RunAction.hh"
#include "Tst16PrimaryGeneratorAction.hh"
#include "Tst16PhysicsList.hh"
#include "Tst16SteppingAction.hh"

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
  runManager->SetUserInitialization(new Tst16DetectorConstruction);
  runManager->SetUserInitialization(new Tst16PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst16RunAction);
  runManager->SetUserAction(new Tst16PrimaryGeneratorAction);
  //runManager->SetUserAction(new Tst16SteppingAction);

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

