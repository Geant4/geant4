// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm5.cc,v 1.4 2000-12-06 18:25:54 maire Exp $
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
#include "Em5VisManager.hh"
#endif

#include "Em5DetectorConstruction.hh"
#include "Em5PhysicsList.hh"
#include "Em5PrimaryGeneratorAction.hh"
#include "Em5RunAction.hh"
#include "Em5EventAction.hh"
#include "Em5SteppingAction.hh"
#include "Em5SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {

  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em5SteppingVerbose);
    
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
  
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
    G4UImanager* UI = G4UImanager::GetUIpointer();  
 
  if (argc==1)   // Define UI terminal for interactive mode  
    { 
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

