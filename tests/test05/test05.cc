// $Id: test05.cc,v 1.3 2000-05-12 14:45:37 gcosmo Exp $
// ------------------------------------------------------------

#include "Tst05DetectorConstruction.hh"
#include "Tst05RunAction.hh"
#include "Tst05PrimaryGeneratorAction.hh"
#include "Tst05PhysicsList.hh"

#ifdef G4VIS_USE
  #include "Tst05VisManager.hh"
#endif

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#if defined(G4UI_USE_XM)
  #include "G4UIXm.hh"
#elif defined(G4UI_USE_XAW)
  #include "G4UIXaw.hh"
#else
  #include "G4UIterminal.hh"
#endif

int main(int argc, char** argv) {

  // Set the default random engine to RanecuEngine
  RanecuEngine defaultEngine;
  HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new Tst05DetectorConstruction);
  runManager->SetUserInitialization(new Tst05PhysicsList);

  // Visualization manager
  #if defined(G4VIS_USE)
    G4VisManager* visManager = new Tst05VisManager;
    visManager->Initialize();
  #endif

  // UserAction classes
  runManager->SetUserAction(new Tst05RunAction);
  runManager->SetUserAction(new Tst05PrimaryGeneratorAction);

  // User interactions
  // Define (G)UI GAG for interactive mode
  if(argc==1)
  {
    // G4UIterminal is a (dumb) terminal
    #if defined(G4UI_USE_XM)
      G4UIsession * session = new G4UIXm(argc,argv);
    #elif defined(G4UI_USE_XAW)
      G4UIsession * session = new G4UIXaw(argc,argv);
    #else
      G4UIsession * session = new G4UIterminal;
    #endif
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

  // Job termination
  #if defined(G4VIS_USE)
    delete visManager;
  #endif
  delete runManager;

  return 0;
}
