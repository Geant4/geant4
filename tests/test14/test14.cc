// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test14.cc,v 1.1 1999-05-29 14:12:02 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TestTst14 
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

#include "Tst14DetectorConstruction.hh"
#include "Tst14PhysicsList.hh"
#include "Tst14PrimaryGeneratorAction.hh"
#include "Tst14RunAction.hh"
#include "Tst14EventAction.hh"
#include "Tst14SteppingAction.hh"

#ifdef GNU_GCC
#include <rw/tpordvec.h>
#include "Tst14CalorHit.hh"
template class RWTPtrOrderedVector <Tst14CalorHit>;
template class RWTPtrVector <Tst14CalorHit>;
template class G4Allocator <Tst14CalorHit>;
#endif

int main(int argc,char** argv) {

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Tst14DetectorConstruction* detector;
  detector = new Tst14DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst14PhysicsList(detector));
  
  // set user action classes
  runManager->SetUserAction(new Tst14PrimaryGeneratorAction(detector));
  Tst14RunAction* runaction = new Tst14RunAction;
  runManager->SetUserAction(runaction);

  Tst14EventAction* eventaction = new Tst14EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Tst14SteppingAction* steppingaction = new Tst14SteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);
    
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
    { 
     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute prerunTst14.mac");    
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

