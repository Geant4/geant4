// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
// 
// --------------------------------------------------------------
//      GEANT 4 - TestEm6 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//      This test is created by V.Ivanchenko on the base of 
//      M.Maire test TestEm5 for the test of LowEnergyIonisation
//      class. Particle is penetrate through the number of sensitive 
//      boxes, dimentions of boxes and particle's kinematic 
//      parameters can be defined by UI interface.
// 29-Jul-1999 V.Ivanchenko first variant 
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
  Em6DetectorConstruction* det;
  runManager->SetUserInitialization(det = new Em6DetectorConstruction);
  runManager->SetUserInitialization(new Em6PhysicsList(det));
  runManager->SetUserAction(new Em6PrimaryGeneratorAction(det));
  
  #ifdef G4VIS_USE
   // visualization manager
   G4VisManager* visManager = new Em6VisManager;
   visManager->Initialize();
  #endif
    
  // set user action classes
  Em6RunAction*   RunAct;
  Em6EventAction* EvtAct;
  
  runManager->SetUserAction(RunAct = new Em6RunAction); 
  runManager->SetUserAction(EvtAct = new Em6EventAction(RunAct));
  //  runManager->SetUserAction(new Em6TrackingAction(RunAct));
  runManager->SetUserAction(new Em6SteppingAction(det,EvtAct,RunAct));
  runManager->SetUserAction(new Em6PrimaryGeneratorAction(det));
  
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







