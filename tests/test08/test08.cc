// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test08.cc,v 1.1 1999-01-08 16:35:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN03
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
#include "T08DetectorConstruction.hh"
#include "T08PhysicsList.hh"
#include "T08RunAction.hh"
#include "T08PrimaryGeneratorAction.hh"
#include "T08EventAction.hh"
#include "T08SteppingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "T08VisManager.hh"
#endif

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  T08DetectorConstruction* T08detector = new T08DetectorConstruction;
  runManager->SetUserInitialization(T08detector);
  runManager->SetUserInitialization(new T08PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new T08VisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new T08RunAction);
  runManager->SetUserAction(new T08PrimaryGeneratorAction(T08detector));
  
  T08EventAction* eventAction = new T08EventAction;
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new T08SteppingAction(T08detector,eventAction));
    
  // User interactions
    G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
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


