// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: hTest.cc,v 1.4 2001-03-12 09:17:14 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TesthTest 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//     
//   
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "hTestVisManager.hh"
#endif

#include "hTestDetectorConstruction.hh"
#include "hTestPhysicsList.hh"
//#include "Tst11PhysicsList.hh"
#include "hTestPrimaryGeneratorAction.hh"
#include "hTestRunAction.hh"
#include "hTestEventAction.hh"
#include "hTestSteppingAction.hh"
#include "hTestSteppingVerbose.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new hTestSteppingVerbose);
    
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  hTestDetectorConstruction* detector;
  detector = new hTestDetectorConstruction;
  runManager->SetUserInitialization(detector);
   runManager->SetUserInitialization(new hTestPhysicsList(detector));
  //runManager->SetUserInitialization(new Tst11PhysicsList());
  /*  
#ifdef G4VIS_USE
  G4cout << "VisManager will be inicialized" << G4endl;
  // visualization manager
  G4VisManager* visManager = new hTestVisManager;
  visManager->Initialize();
#endif 
*/
  G4cout << "User actions will be inicialized" << G4endl;
 
  // set user action classes
  runManager->SetUserAction(new hTestPrimaryGeneratorAction(detector));
  hTestRunAction* runaction = new hTestRunAction;
  runManager->SetUserAction(runaction);

  hTestEventAction* eventaction = new hTestEventAction(runaction);
  runManager->SetUserAction(eventaction);

  hTestSteppingAction* steppingaction = new hTestSteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
    { 
     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute init.mac");    
     session->SessionStart();
     delete session;
    }
  else           // Batch mode
    { 
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }
    
  // job termination
  /*
#ifdef G4VIS_USE
  delete visManager;
#endif
*/
  G4cout << "runManager will be deleted" << G4endl;  
  delete runManager;
  G4cout << "runManager is deleted" << G4endl;  

  return 0;

}

