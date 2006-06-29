//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: testNTST.cc,v 1.3 2006-06-29 18:25:19 gunter Exp $
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

