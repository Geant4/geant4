// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleN02.cc,v 1.1 1999-01-07 16:05:43 gunter Exp $
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

#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02RunAction.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "ExN02VisManager.hh"
#endif

#include "G4ios.hh"


int
main(int argc,char** argv) {

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  ExN02DetectorConstruction* ExN02detector = new ExN02DetectorConstruction;
  runManager->SetUserInitialization(ExN02detector);
  runManager->SetUserInitialization(new ExN02PhysicsList);
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new ExN02VisManager;
  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new ExN02RunAction);
  runManager->SetUserAction(new ExN02PrimaryGeneratorAction(ExN02detector));
  
  ExN02EventAction* eventAction = new ExN02EventAction;
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new ExN02SteppingAction(ExN02detector,eventAction));
    
  // User interactions
    G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
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


