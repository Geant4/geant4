// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm2.cc,v 1.1 1999-02-19 10:28:37 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#include "Em2VisManager.hh"
#endif

#include "Em2DetectorConstruction.hh"
#include "Em2PhysicsList.hh"
#include "Em2PrimaryGeneratorAction.hh"
#include "Em2RunAction.hh"
#include "Em2EventAction.hh"
#include "Em2TrackingAction.hh"
#include "Em2SteppingAction.hh"

int main(int argc,char** argv) {

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em2DetectorConstruction* detector = new Em2DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em2PhysicsList);
  
  Em2PrimaryGeneratorAction* primary = new Em2PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
    
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em2VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  Em2RunAction* RunAct = new Em2RunAction(detector,primary);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Em2EventAction   (RunAct));
  runManager->SetUserAction(new Em2TrackingAction(RunAct));
  runManager->SetUserAction(new Em2SteppingAction(detector,RunAct));
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode.
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

