// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testAssemblyVolume.cc,v 1.2 2001-02-01 21:24:28 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT4 - testAssemblyVolume 
//
//      For information related to this code contact:
//      CERN, IT Division, API Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "TstVADetectorConstruction.hh"
#include "TstVARunAction.hh"
#include "TstVAEventAction.hh"
#include "TstVAPrimaryGeneratorAction.hh"
#include "TstVASteppingAction.hh"
#include "TstVAPhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#ifdef G4VIS_USE
#include "TstVAVisManager.hh"
#endif

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new TstVADetectorConstruction);
  runManager->SetUserInitialization(new TstVAPhysicsList);

  #ifdef G4VIS_USE
    // Visualization, if you choose to have it!
    G4VisManager* visManager = new TstVAVisManager;
    visManager->Initialize();
  #endif

  // UserAction classes
  runManager->SetUserAction(new TstVARunAction);
  runManager->SetUserAction(new TstVAPrimaryGeneratorAction);
  runManager->SetUserAction(new TstVASteppingAction);
  runManager->SetUserAction(new TstVAEventAction);

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

