// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testDrawVox.cc,v 1.1 1999-07-28 17:57:21 graignac Exp $
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

#include "TstDrawVox01DetectorConstruction.hh"
#include "TstDrawVox01RunAction.hh"
#include "TstDrawVox01PrimaryGeneratorAction.hh"
#include "TstDrawVox01SteppingAction.hh"
#include "TstDrawVox01PhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "TstDrawVox01VisManager.hh"
#endif

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new TstDrawVox01DetectorConstruction);
  runManager->SetUserInitialization(new TstDrawVox01PhysicsList);

  #ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    G4VisManager* visManager = new TstDrawVox01VisManager;
    visManager->Initialize();
  #endif

  // UserAction classes
  runManager->SetUserAction(new TstDrawVox01RunAction);
  runManager->SetUserAction(new TstDrawVox01PrimaryGeneratorAction);
  runManager->SetUserAction(new TstDrawVox01SteppingAction);

  // User interactions
  // Define (G)UI for interactive mode
 
  // User interactions
    G4UImanager * UI = G4UImanager::GetUIpointer();
  
  if(argc==1)
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerun.g4mac"); 
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
  
  #ifdef G4VIS_USE
    delete visManager;
  #endif
  delete runManager;
  return 0;
}

