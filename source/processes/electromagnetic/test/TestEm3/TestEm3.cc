// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm3.cc,v 1.2 1999-01-26 16:27:34 maire Exp $
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
#include "Em3VisManager.hh"
#endif

#include "Em3DetectorConstruction.hh"
#include "Em3PhysicsList.hh"
#include "Em3PrimaryGeneratorAction.hh"
#include "Em3RunAction.hh"
#include "Em3EventAction.hh"
#include "Em3SteppingAction.hh"

#ifdef GNU_GCC
#include <rw/tpordvec.h>
#include "Em3CalorHit.hh"
template class RWTPtrOrderedVector <Em3CalorHit>;
template class RWTPtrVector <Em3CalorHit>;
template class G4Allocator  <Em3CalorHit>;
#endif

int main(int argc,char** argv) {

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em3DetectorConstruction* detector = new Em3DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em3PhysicsList);
  
  runManager->SetUserAction(new Em3PrimaryGeneratorAction(detector));
    
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em3VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  Em3RunAction* RunAct = new Em3RunAction;
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Em3EventAction(RunAct));
  runManager->SetUserAction(new Em3SteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
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

