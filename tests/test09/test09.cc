// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test09.cc,v 1.1 1999-01-08 16:35:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst09DetectorConstruction.hh"
#include "Tst09RunAction.hh"
#include "Tst09PrimaryGeneratorAction.hh"
#include "Tst09PhysicsList.hh"
#include "Tst09SteppingAction.hh"
#include "Tst09TrackingAction.hh"

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
  runManager->SetUserInitialization(new Tst09DetectorConstruction);
  runManager->SetUserInitialization(new Tst09PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst09RunAction);
  runManager->SetUserAction(new Tst09PrimaryGeneratorAction);
  //runManager->SetUserAction(new Tst09SteppingAction);
  runManager->SetUserAction(new Tst09TrackingAction);

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

