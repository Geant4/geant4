// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTest ----------------------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "hTestVisManager.hh"
#endif

#include "hTestDetectorConstruction.hh"
#include "hTestPhysicsList.hh"
#include "hTestPrimaryGeneratorAction.hh"
#include "hTestRunAction.hh"
#include "hTestEventAction.hh"
#include "hTestTrackingAction.hh"
#include "hTestSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) {

  //choose the Random engine
  //  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager();

  // set mandatory initialization classes
  hTestDetectorConstruction* det = new hTestDetectorConstruction();
  runManager->SetUserInitialization(det);
  det->SetVerbose(2);
  G4cout << "Detector Construction is defined" << G4endl;
  
  runManager->SetUserInitialization(new hTestPhysicsList(det));
  G4cout << "Physics List is defined" << G4endl;

#ifdef G4VIS_USE
  G4cout << "VisManager will be inicialized" << G4endl;
  // visualization manager
  G4VisManager* visManager = new hTestVisManager;
  visManager->Initialize();
#endif 

  G4cout << "User actions will be initialized" << G4endl;
 
  // set user action classes
  runManager->SetUserAction(new hTestPrimaryGeneratorAction(det));
  hTestRunAction* run = new hTestRunAction(det);
  runManager->SetUserAction(run);
  G4cout << "RunAction is defined" << G4endl;

  hTestEventAction* event = new hTestEventAction(run,det);
  runManager->SetUserAction(event);
  G4cout << "EventAction is defined" << G4endl;

  runManager->SetUserAction(new hTestTrackingAction(run));
  G4cout << "TrackingAction is defined" << G4endl;
  runManager->SetUserAction(new hTestSteppingAction(run,event,det));
  G4cout << "SteppingAction is defined" << G4endl;
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  
     
  if (argc==1)   // Define UI terminal for interactive mode  
    { 
     // Initialize G4 kernel
     runManager->Initialize();

     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute init.mac");    
     session->SessionStart();
     delete session;
    }
  else  // Batch mode with 1 or more files no beamOn command in the 1st
    { 
     G4cout << "UI interface is started" << G4endl;
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);

     // Initialize G4 kernel
     G4cout << "Start initialisation" << G4endl;
     runManager->Initialize();

     // 
     G4cout << "Start event loop" << G4endl;
     runManager->BeamOn(det->GetNumberOfEvents());

     // next file
     if(argc>1) UI->ApplyCommand(command+argv[2]);
    }
    
  // job termination
  
#ifdef G4VIS_USE
  delete visManager;
#endif

  G4cout << "runManager will be deleted" << G4endl;  
  delete runManager;
  G4cout << "runManager is deleted" << G4endl;  

  return 0;
}

