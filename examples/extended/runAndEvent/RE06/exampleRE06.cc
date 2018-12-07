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
/// \file RE06/exampleRE06.cc
/// \brief Main program of the RE06 example
//
//
// --------------------------------------------------------------
//      GEANT 4 - example RE06
// --------------------------------------------------------------

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4ParallelWorldPhysics.hh"

#include "RE06DetectorConstruction.hh"
#include "RE06ParallelWorld.hh"
#include "RE06PrimaryGeneratorAction.hh"
#include "RE06RunAction.hh"
#include "RE06SteppingVerbose.hh"
#include "RE06ActionInitialization.hh"
#include "RE06WorkerInitialization.hh"

int main(int argc,char** argv)
{
 // Instantiate G4UIExecutive if there are no arguments (interactive mode)
 G4UIExecutive* ui = nullptr;
 if ( argc == 1 ) {
   ui = new G4UIExecutive(argc, argv);
 }

 // Construct the stepping verbose class
 RE06SteppingVerbose* verbosity = new RE06SteppingVerbose;
 G4VSteppingVerbose::SetInstance(verbosity);

 // Construct the run manager
 //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4);
  runManager->SetUserInitialization(new RE06WorkerInitialization);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

 // Set mandatory initialization classes
 //
 G4String parallelWorldName = "ParallelScoringWorld";
 G4VUserDetectorConstruction* detector = new RE06DetectorConstruction;
 detector->RegisterParallelWorld(new RE06ParallelWorld(parallelWorldName));
 runManager->SetUserInitialization(detector);
 //
 G4VModularPhysicsList* physics = new FTFP_BERT;
 physics->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));
 runManager->SetUserInitialization(physics);

 // Set user action classes
 //
 runManager->SetUserInitialization(new RE06ActionInitialization);

 // Visualization manager
 //
 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();

 // Initialize G4 kernel
 //
 runManager->Initialize();

 // Get the pointer to the User Interface manager
 //
 G4UImanager* UImanager = G4UImanager::GetUIpointer();

 if (ui)   // Define UI session for interactive mode
 {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
 }
 else           // Batch mode
 {
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   UImanager->ApplyCommand(command+fileName);
 }

  // Job termination
  // Free the store:
  // user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not
  // be deleted in the main() program !

 delete visManager;
 delete runManager;

 return 0;
}
