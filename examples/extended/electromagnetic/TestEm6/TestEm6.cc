// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm6.cc,v 1.6 2000-01-28 04:25:10 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TestEm6 
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
#include "Em6VisManager.hh"
#endif

#include "Em6DetectorConstruction.hh"
#include "Em6PhysicsList.hh"
#include "Em6PrimaryGeneratorAction.hh"
#include "Em6RunAction.hh"
#include "Em6EventAction.hh"
#include "Em6SteppingAction.hh"
#include "Em6SteppingVerbose.hh"

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em6SteppingVerbose);
    
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em6DetectorConstruction* detector;
  detector = new Em6DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em6PhysicsList(detector));
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em6VisManager;
  visManager->Initialize();
#endif 
 
  // set user action classes
  runManager->SetUserAction(new Em6PrimaryGeneratorAction(detector));
  Em6RunAction* runaction = new Em6RunAction;
  runManager->SetUserAction(runaction);

  Em6EventAction* eventaction = new Em6EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Em6SteppingAction* steppingaction = new Em6SteppingAction(detector,
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

