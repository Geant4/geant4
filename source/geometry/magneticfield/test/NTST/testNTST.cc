// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testNTST.cc,v 1.1 2003-11-27 09:02:46 japost Exp $
// 
// ----------------------------------------------------------------
//      GEANT 4 - exampleNTST : BaBar SVT standalone simulation
//
//      For information related to this code contact: Bill Lockman
//                                       lockman@slac.stanford.edu
// ----------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
// #include "g4std/stdio.h"

#ifdef G4VIS_USE
#include "NTSTVisManager.hh"
#endif

#include "NTSTFileRead.hh"
#include "NTSTDetectorConstruction.hh"
#include "NTSTPhysicsList.hh"
#include "NTSTPrimaryGeneratorAction.hh"
#include "NTSTRunAction.hh"
#include "NTSTEventAction.hh"
#include "NTSTSteppingAction.hh"

#ifdef TRAP_ON_NAN_ENABLE
#include <stdio.h>
#include <fenv.h>
#endif

void trap_on_NaN ()
{
#ifdef TRAP_ON_NAN_ENABLE
  // fenv_t  ALL_Except(FE_ALL_EXCEPT);  // fails
  fenv_t fenv_state;

  fegetenv(&fenv_state); 
  fenv_state |= FE_ALL_EXCEPT;
  // fesetenv (FE_NOMASK_ENV);
  // G4cout << "Trying to set trap (signal) for NaN: using FE_NOMASK_ENV" << G4endl;

  // fesetenv (FE_DIVBYZERO|FE_INVALID);
  // G4cout << "Trying to set trap (signal) for NaN: using FE_DIVBYZERO and FE_INVALID" << G4endl;

  // fesetenv ((fenv_t *)FE_ALL_EXCEPT); 
  fesetenv (ALL_Except); 
  // fesetenv (&ALL_Except); 
  G4cout << "Trying to set trap (signal) for NaN: using FE_ALL_EXCEPT" << G4endl;
  

  // fesetenv (FE_DIVBYZERO); 
  // G4cout << "Trying to set trap (signal) for NaN: using FE_DIVBYZERO" << G4endl;

  // G4cout << "Setting trap for NaN using fpu_" << G4endl;
  // fpu_();
#endif
}

#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

int main(int argc,char** argv) {

    trap_on_NaN();

        // Construct the default run manager
    G4RunManager * runManager = new G4RunManager;
    
        // set mandatory initialization classes
    NTSTDetectorConstruction* detector = new NTSTDetectorConstruction();
    NTSTPhysicsList *physicsList = new NTSTPhysicsList();
  
        // detector construction object
  runManager->SetUserInitialization(detector);
  
  runManager->SetUserInitialization(physicsList);
 
  // G4PropagatorInField* propagator=
  //   G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
  // propagator->SetVerboseLevel(3); 

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new NTSTVisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new NTSTPrimaryGeneratorAction);
  runManager->SetUserAction(new NTSTRunAction);
  runManager->SetUserAction(new NTSTEventAction);
  runManager->SetUserAction(new NTSTSteppingAction);
    
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode  
    { 
      G4UIsession * session = new G4UIterminal(new G4UItcsh);
      //      UI->ApplyCommand("/control/execute prerunNTST.mac");    
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

