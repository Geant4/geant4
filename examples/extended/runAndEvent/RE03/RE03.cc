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
/// \file runAndEvent/RE03/RE03.cc
/// \brief Main program of the runAndEvent/RE03 example
//
//
//

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"

#include "RE03DetectorConstruction.hh"
#include "FTFP_BERT.hh" 
#include "RE03ActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//====================================================================
// Un-comment this line for user defined score writer
//    #include "RE03UserScoreWriter.hh"
//====================================================================

int main(int argc,char** argv)
{
 // Instantiate G4UIExecutive if there are no arguments (interactive mode)
 G4UIExecutive* ui = nullptr;
 if ( argc == 1 ) {
   ui = new G4UIExecutive(argc, argv);
 }

 // Construct the run manager
 //
#ifdef G4MULTITHREADED
 G4MTRunManager * runManager = new G4MTRunManager;
 runManager->SetNumberOfThreads(4);
#else
 G4RunManager * runManager = new G4RunManager;
#endif

 // Activate UI-command base scorer
 G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
 scManager->SetVerboseLevel(1);

//====================================================================
// Un-comment this line for user defined score writer
//    scManager->SetScoreWriter(new RE03UserScoreWriter());
//====================================================================

 // Set mandatory initialization classes
 //
 G4VUserDetectorConstruction* detector = new RE03DetectorConstruction;
 runManager->SetUserInitialization(detector);
 //
 G4VUserPhysicsList* physics = new FTFP_BERT;
 runManager->SetUserInitialization(physics);
    
 // Set user action classes through Worker Initialization
 //
 RE03ActionInitialization* actions = new RE03ActionInitialization;
 runManager->SetUserInitialization(actions);
  
 // Visualization manager
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
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

 delete visManager;
 delete runManager;

 return 0;
}
