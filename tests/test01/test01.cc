// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test01.cc,v 1.2 1999-11-26 09:46:58 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT4 - test01 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "Tst01DetectorConstruction.hh"
#include "Tst01RunAction.hh"
#include "Tst01PrimaryGeneratorAction.hh"
#include "Tst01SteppingAction.hh"
#include "Tst01PhysicsList.hh"

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
  runManager->SetUserInitialization(new Tst01DetectorConstruction);
  runManager->SetUserInitialization(new Tst01PhysicsList);

  // UserAction classes
  runManager->SetUserAction(new Tst01RunAction);
  runManager->SetUserAction(new Tst01PrimaryGeneratorAction);
  runManager->SetUserAction(new Tst01SteppingAction);

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

