// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm0.cc,v 1.3 1999-05-10 16:44:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - example Em0 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

///#define Em0NoOptimize 1

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4UIterminal.hh"

#include "Em0DetectorConstruction.hh"
#include "Em0PhysicsList.hh"
#include "Em0PrimaryGeneratorAction.hh"

#ifdef Em0NoOptimize
#include "Em0RunAction.hh"
#include "Em0EventAction.hh"
#include "Em0TrackingAction.hh"
#include "Em0SteppingAction.hh"
#endif

int main(int argc,char** argv) {

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em0DetectorConstruction* DetConst;
  runManager->SetUserInitialization(DetConst = new Em0DetectorConstruction);
  runManager->SetUserInitialization(new Em0PhysicsList(DetConst));
  runManager->SetUserAction(new Em0PrimaryGeneratorAction);
  
#ifdef Em0NoOptimize
  // set user action classes
  Em0RunAction*   RunAct;
  Em0EventAction* EvtAct;
  
  runManager->SetUserAction(RunAct = new Em0RunAction); 
  runManager->SetUserAction(EvtAct = new Em0EventAction);
  runManager->SetUserAction(new Em0TrackingAction(RunAct));
  runManager->SetUserAction(new Em0SteppingAction(EvtAct));
#endif
  
  //Initialize G4 kernel
  //runManager->Initialize();
    
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
  delete runManager;

  return 0;
}

