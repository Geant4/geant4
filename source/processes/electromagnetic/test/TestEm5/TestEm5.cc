// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm5.cc,v 1.1 1999-01-08 16:33:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - TestEm5 
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

#ifdef G4VIS_USE
#include "Em5VisManager.hh"
#endif

#include "Em5DetectorConstruction.hh"
#include "Em5PhysicsList.hh"
#include "Em5PrimaryGeneratorAction.hh"
#include "Em5RunAction.hh"
#include "Em5EventAction.hh"
#include "Em5SteppingAction.hh"

// Solve templates.
// Code in g4templates.hh controlled by the macro G4_SOLVE_TEMPLATES
#ifdef G4_SOLVE_TEMPLATES
#ifdef G4VIS_USE
#  define G4_SOLVE_VIS_TEMPLATES
#endif
#endif
#include "g4templates.hh"

#ifdef GNU_GCC
#include <rw/tpordvec.h>
#include "Em5CalorHit.hh"
template class RWTPtrOrderedVector <Em5CalorHit>;
template class RWTPtrVector <Em5CalorHit>;
template class G4Allocator <Em5CalorHit>;
#endif

int main(int argc,char** argv) {

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em5DetectorConstruction* detector;
  detector = new Em5DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em5PhysicsList(detector));
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em5VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new Em5PrimaryGeneratorAction(detector));
  Em5RunAction* runaction = new Em5RunAction;
  runManager->SetUserAction(runaction);

  Em5EventAction* eventaction = new Em5EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Em5SteppingAction* steppingaction = new Em5SteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);
    
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
    { 
     G4UIsession * session = new G4UIterminal;
     UI->ApplyCommand("/control/execute prerunEm5.mac");    
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

