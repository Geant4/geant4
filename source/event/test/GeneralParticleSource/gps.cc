// 
// This code has been created by C Ferguson, based on
// example GEANT4 code, for a University of Southampton/
// DERA/ESA/CERN collaboration. Further information
// can be obtained from C Ferguson at the University of 
// Southampton.
// 

#include "ExampleDetectorConstruction.hh"
#include "ExamplePhysicsList.hh"
#include "ExamplePrimaryGeneratorAction.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
//#include "ExampleVisManager.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  int seed, count;

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  ExampleDetectorConstruction *Mydetector = new ExampleDetectorConstruction;
  runManager->SetUserInitialization(Mydetector);
  runManager->SetUserInitialization(new ExamplePhysicsList);

  // Visualization, if you choose to have it!
  //  G4VisManager* visManager = new ExampleVisManager;
  //visManager->Initialize();
   
  // UserAction classes
  runManager->SetUserAction(new ExamplePrimaryGeneratorAction);
  
  // User interactions
  G4UImanager * UI = G4UImanager::GetUIpointer();  
  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/event/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");

  if(argc==1)
    // Define (G)UI terminal for interactive mode  
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = new G4UIterminal;
      UI->ApplyCommand("/control/execute prerun.g4mac");    
      session->SessionStart();
      delete session;
    }
  else {
    // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  //delete visManager;
  delete runManager;

  return 0;
}








