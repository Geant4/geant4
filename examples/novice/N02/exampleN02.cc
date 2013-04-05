//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN02DetectorConstruction.hh"
#include "ExN02PhysicsList.hh"
#include "ExN02PrimaryGeneratorAction.hh"
#include "ExN02RunAction.hh"
#include "ExN02EventAction.hh"
#include "ExN02SteppingAction.hh"
#include "ExN02SteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"

#ifdef G4MULTITHREADED
#include "ExN02WorkerInitialization.hh"
#include "G4MTRunManager.hh"
#endif

#include "FTFP_BERT.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
int main(int argc,char** argv)
{
    //TODO: To remove
    CLHEP::RanluxEngine defaultEngine( 1234567, 4 );
    G4Random::setTheEngine( &defaultEngine );
    //G4int seed = time( NULL );
    //  CLHEP::HepRandom::setTheSeed( seed );
    G4Random::setTheSeed( 1220515164 );
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new ExN02SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  
  // Run manager
  //
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    runManager->SetNumberThreads( 2 );
    //Common to all threads
    ExN02DetectorConstruction* detector = new ExN02DetectorConstruction;
    runManager->SetUserInitialization(detector);
    ExN02WorkerInitialization* workerInit = new ExN02WorkerInitialization;
    runManager->SetUserInitialization( workerInit );

    
    //The following three lines are specific of this application
    workerInit->SetDetectorConstruction(detector);
    G4String fileName = argv[1];
    workerInit->SetMacroFileName(fileName);
    
    //Needed for initialization...
    G4VUserPhysicsList* physics = new ExN02PhysicsList;
    //G4VUserPhysicsList* physics = new FTFP_BERT;
    runManager->SetUserInitialization(physics);

    // User Action classes., These should not be needed
    //Problem: this action defines the /gun commands that, if used in macroFile used below
    //will cause the job to stop.
    //Why we need to execute the macroFile? Because it may contains commands that modify the
    //setup (e.g. geometry related).
    
    //TODO: Implement a correct way that the "master" thread can skip commands that are not found
    //G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
    //srunManager->SetUserAction(gen_action);
    //
    //G4UserRunAction* run_action = new ExN02RunAction;
    //runManager->SetUserAction(run_action);
    //
    //G4UserEventAction* event_action = new ExN02EventAction;
    //runManager->SetUserAction(event_action);
    //
    //G4UserSteppingAction* stepping_action = new ExN02SteppingAction;
    //runManager->SetUserAction(stepping_action);
    
    
    // Initialize G4 kernel
    //
    runManager->Initialize();
    
    // Get the pointer to the User Interface manager
    //
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    
    if (argc!=1)   // batch mode
    {
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command+fileName);
    }
    else           // interactive mode : define UI session
    {
        //Not supported yet!!!!
        return 1;
    }
    
    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
#else
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  ExN02DetectorConstruction* detector = new ExN02DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new ExN02PhysicsList;
  runManager->SetUserInitialization(physics);
   
  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new ExN02PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  G4UserRunAction* run_action = new ExN02RunAction;
  runManager->SetUserAction(run_action);
  //
  G4UserEventAction* event_action = new ExN02EventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new ExN02SteppingAction;
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    
     
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
#endif //G4MULTITHREADED
  delete runManager;
  delete verbosity;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

