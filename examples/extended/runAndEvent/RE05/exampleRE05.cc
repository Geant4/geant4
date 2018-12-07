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
/// \file RE05/exampleRE05.cc
/// \brief Main program of the RE05 example
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "RE05WorkerInitialization.hh"
#else
#include "G4RunManager.hh"
#include "RE05SteppingVerbose.hh"
#endif

#include "G4UImanager.hh"

#include "RE05DetectorConstruction.hh"
#include "RE05CalorimeterParallelWorld.hh"
#include "QBBC.hh"
#include "G4ParallelWorldPhysics.hh"
#include "RE05ActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  G4int number_of_threads = 4; // default number of threads
  runManager->SetNumberOfThreads(number_of_threads);
  runManager->SetUserInitialization(new RE05WorkerInitialization);
#else
  G4VSteppingVerbose* verbosity = new RE05SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);
  G4RunManager* runManager = new G4RunManager;
#endif

  G4String parallelWorldName = "ReadoutWorld";
  // User Initialization classes (mandatory)
  //
  G4VUserDetectorConstruction* detector = new RE05DetectorConstruction();
  detector->RegisterParallelWorld
       (new RE05CalorimeterParallelWorld(parallelWorldName));
  runManager->SetUserInitialization(detector);
  //
  G4VModularPhysicsList* physicsList = new QBBC;
  physicsList
   ->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));
  runManager->SetUserInitialization(physicsList);
  //
  G4VUserActionInitialization* actions = new RE05ActionInitialization;
  runManager->SetUserInitialization(actions);

  runManager->Initialize();

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  //get the pointer to the User Interface manager   
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (!ui)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define UI session
    {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
    }

  delete visManager;

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
