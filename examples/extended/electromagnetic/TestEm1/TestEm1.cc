// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TestEm1.cc,v 1.4 2000-12-06 13:21:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

#define Em1NoOptimize 1

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "Em1DetectorConstruction.hh"
#include "Em1PhysicsList.hh"
#include "Em1PrimaryGeneratorAction.hh"
#include "Em1SteppingVerbose.hh"

#ifdef Em1NoOptimize
 #include "Em1RunAction.hh"
 #include "Em1EventAction.hh"
 #include "Em1TrackingAction.hh"
 #include "Em1SteppingAction.hh"

 #ifdef G4VIS_USE
  #include "Em1VisManager.hh"
 #endif
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
 
int main(int argc,char** argv) {
 
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em1SteppingVerbose);
    
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em1DetectorConstruction* det;
  runManager->SetUserInitialization(det = new Em1DetectorConstruction);
  runManager->SetUserInitialization(new Em1PhysicsList(det));
  runManager->SetUserAction(new Em1PrimaryGeneratorAction);
  
#ifdef Em1NoOptimize  
  #ifdef G4VIS_USE
   // visualization manager
   G4VisManager* visManager = new Em1VisManager;
   visManager->Initialize();
  #endif
    
  // set user action classes
  Em1RunAction*   RunAct;
  Em1EventAction* EvtAct;
  
  runManager->SetUserAction(RunAct = new Em1RunAction); 
  runManager->SetUserAction(EvtAct = new Em1EventAction(RunAct));
  runManager->SetUserAction(new Em1TrackingAction(RunAct));
  runManager->SetUserAction(new Em1SteppingAction(RunAct,EvtAct));
#endif
  
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
#ifdef Em1NoOptimize   
 #ifdef G4VIS_USE
  delete visManager;
 #endif
#endif
 
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
