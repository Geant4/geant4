// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm3.cc,v 1.4 2000-12-06 16:43:38 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "Em3VisManager.hh"
#endif

#include "Em3DetectorConstruction.hh"
#include "Em3PhysicsList.hh"
#include "Em3PrimaryGeneratorAction.hh"
#include "Em3RunAction.hh"
#include "Em3EventAction.hh"
#include "Em3SteppingAction.hh"
#include "Em3SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em3SteppingVerbose);
  
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
  Em3RunAction* RunAct = new Em3RunAction(detector);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Em3EventAction(RunAct,detector));
  runManager->SetUserAction(new Em3SteppingAction);
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif                 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

