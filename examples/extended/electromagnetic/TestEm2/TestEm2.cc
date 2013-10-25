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
/// \file electromagnetic/TestEm2/TestEm2.cc
/// \brief Main program of the electromagnetic/TestEm2 example
//
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "ActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);


 // Number of threads is defined via 3nd argument
  G4String nn = "";

#ifdef G4MULTITHREADED  
  if (argc==3) { nn = argv[2]; }

  if("" == nn) {
    // Number of threads is defined via environment variable
    char* path = getenv("G4NUMBEROFTHREADS");
    if(path) { nn = G4String(path); }
  }
#endif

  G4RunManager * runManager = 0;
  G4MTRunManager * runManagerMT = 0;

  // MT mode
  if("" != nn) {
    runManagerMT = new G4MTRunManager();
    G4int N = 0;
    std::istringstream is(nn);
    is >> N;
    if(N < 1) { N = 1; }
    runManagerMT->SetNumberOfThreads(N);

    // set mandatory initialization classes
    DetectorConstruction* detector = new DetectorConstruction();
    runManagerMT->SetUserInitialization(detector);
    runManagerMT->SetUserInitialization(new PhysicsList());
    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);

    // set user actions
    runManagerMT->SetUserInitialization(new ActionInitialization(detector,
                                                                 primary));

    G4cout << "##### TestEm2 started for " << N << " threads" 
           << " #####" << G4endl;

    // sequential mode
  } else {
    runManager = new G4RunManager();

    // set mandatory initialization classes
    DetectorConstruction* detector = new DetectorConstruction();
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList());  
    PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
    runManager->SetUserAction(primary);
    
    // set user action classes
    RunAction* RunAct = new RunAction(detector,primary);
    runManager->SetUserAction(RunAct);
    runManager->SetUserAction(new EventAction   (RunAct));
    runManager->SetUserAction(new TrackingAction(RunAct));
    runManager->SetUserAction(new SteppingAction(detector,RunAct)); 

    G4cout << "##### test40 started in sequential mode" 
           << " #####" << G4endl;

  }

  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
    
  else  // define visualization and UI terminal for interactive mode
    { 
#ifdef G4VIS_USE
     G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
#endif
                 
     
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
#endif
      
#ifdef G4VIS_USE
      delete visManager;
#endif            
    }

  // job termination
  //
  delete runManager;
  delete runManagerMT;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
