// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * main program class                                               *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "LISAVisManager.hh"
#endif

#include "LISADetectorConstruction.hh"
#include "LISAPhysicsList.hh"
#include "LISAPrimaryGeneratorAction.hh"
#include "LISARunAction.hh"
#include "LISAEventAction.hh"
#include "LISASteppingAction.hh"
#include "LISAStackingAction.hh"

#include <vector>


int main(int argc,char** argv) {
  
  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes
  LISADetectorConstruction* detector = new LISADetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new LISAPhysicsList);
  LISAPrimaryGeneratorAction* generatorAction = new LISAPrimaryGeneratorAction;
  runManager->SetUserAction(generatorAction);
  
  // set user action classes
  LISASteppingAction* steppingAction = new LISASteppingAction;
  runManager->SetUserAction(steppingAction);
  runManager->SetUserAction(new LISAStackingAction);
  runManager->SetUserAction(
     new LISAEventAction(generatorAction,steppingAction));
  runManager->SetUserAction(new LISARunAction);

  
  G4UIsession* session=0;
  if (argc==1) {  // Define UI session for interactive mode.
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
    session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif
#endif
  }
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new LISAVisManager;
  visManager->Initialize();
#endif
  

  //Initialize G4 kernel
  runManager->Initialize();
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  
  // Define UI session for interactive mode.
  if(session) {  
    // G4UIterminal is a (dumb) terminal.
    UI->ApplyCommand("/control/execute init.mac");    
    session->SessionStart();
    delete session;
  }
  // Batch mode
  else { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }
  
  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

