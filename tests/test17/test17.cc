// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test17.cc,v 1.1 1999-11-30 18:01:49 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - test17 based on TestEm6
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

#include "Tst17DetectorConstruction.hh"
#include "Tst17PhysicsList.hh"
#include "Tst17PrimaryGeneratorAction.hh"
#include "Tst17RunAction.hh"
#include "Tst17EventAction.hh"
#include "Tst17SteppingAction.hh"

#ifdef G4VIS_USE
#include "Tst17VisManager.hh"
#endif

#ifdef GNU_GCC
#include "g4rw/tpordvec.h"
#include "Tst17CalorHit.hh"
template class G4RWTPtrOrderedVector <Tst17CalorHit>;
template class G4RWTPtrVector <Tst17CalorHit>;
template class G4Allocator <Tst17CalorHit>;
#endif

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Tst17DetectorConstruction* detector;
  detector = new Tst17DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Tst17PhysicsList(detector));
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Tst17VisManager;
  visManager->Initialize();
#endif

  // set user action classes
  runManager->SetUserAction(new Tst17PrimaryGeneratorAction(detector));
  Tst17RunAction* runaction = new Tst17RunAction;
  runManager->SetUserAction(runaction);

  Tst17EventAction* eventaction = new Tst17EventAction(runaction);
  runManager->SetUserAction(eventaction);

  Tst17SteppingAction* steppingaction = new Tst17SteppingAction(detector,
                                               eventaction, runaction);
  runManager->SetUserAction(steppingaction);

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

#ifdef G4VIS_USE
  delete visManager;
#endif

  // job termination
  delete runManager;

  return 0;
}

