// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleN03.cc,v 1.2 1999-03-03 14:53:34 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN02 
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
#include "ExN03VisManager.hh"
#endif

#include "ExN03DetectorConstruction.hh"
#include "ExN03PhysicsList.hh"
#include "ExN03PrimaryGeneratorAction.hh"
#include "ExN03RunAction.hh"
#include "ExN03EventAction.hh"
#include "ExN03SteppingAction.hh"


#ifdef GNU_GCC
#include <rw/tpordvec.h>
#include "ExN03CalorHit.hh"
template class RWTPtrOrderedVector <ExN03CalorHit>;
template class RWTPtrVector <ExN03CalorHit>;
template class G4Allocator <ExN03CalorHit>;
#endif

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  ExN03DetectorConstruction* detector = new ExN03DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new ExN03PhysicsList);
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new ExN03VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  runManager->SetUserAction(new ExN03PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new ExN03RunAction);
  runManager->SetUserAction(new ExN03EventAction);
  runManager->SetUserAction(new ExN03SteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = new G4UIterminal;
      UI->ApplyCommand("/control/execute prerunN03.mac");    
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

