// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test17.cc,v 1.4 2000-05-26 15:50:12 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TestTest17 
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
#include "Test17VisManager.hh"
#endif

#include "Test17DetectorConstruction.hh"
#include "Test17PhysicsList.hh"
#include "Test17PrimaryGeneratorAction.hh"
#include "Test17RunAction.hh"
#include "Test17EventAction.hh"
#include "Test17SteppingAction.hh"
#include "Test17SteppingVerbose.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Test17SteppingVerbose);
    
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Test17DetectorConstruction* detector;
  detector = new Test17DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Test17PhysicsList(detector));
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Test17VisManager;
  visManager->Initialize();
#endif 
 
  // set user action classes
  runManager->SetUserAction(new Test17PrimaryGeneratorAction(detector));
  Test17RunAction* runaction = new Test17RunAction;
  runManager->SetUserAction(runaction);

  Test17EventAction* eventaction = new Test17EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Test17SteppingAction* steppingaction = new Test17SteppingAction(detector,
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
#ifdef G4VIS_USE
  delete visManager;
#endif  
  delete runManager;

  return 0;
}

