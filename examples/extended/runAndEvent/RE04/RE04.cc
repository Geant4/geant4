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
/// \file runAndEvent/RE04/RE04.cc
/// \brief Main program of the runAndEvent/RE04 example
//
// $Id: $
//
#include "G4RunManager.hh"
#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "RE04DetectorConstruction.hh"
#include "RE04ParallelWorldConstruction.hh"
#include "RE04PhysicsList.hh"
#include "RE04PrimaryGeneratorAction.hh"
#include "RE04EventAction.hh"
#include "RE04TrackingAction.hh"
#include "RE04SteppingAction.hh"

int main(int argc,char** argv)
{
 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;
 G4ScoringManager::GetScoringManager();

 G4String paraWorldName = "ParallelWorld";

 // Set mandatory initialization classes
 //
 RE04DetectorConstruction* realWorld = new RE04DetectorConstruction;
 RE04ParallelWorldConstruction* parallelWorld
   = new RE04ParallelWorldConstruction(paraWorldName);
 realWorld->RegisterParallelWorld(parallelWorld);
 runManager->SetUserInitialization(realWorld);
 //
 G4VUserPhysicsList* physics = new RE04PhysicsList(paraWorldName);
 runManager->SetUserInitialization(physics);
    
 // Set user action classes
 //
 runManager->SetUserAction(new RE04PrimaryGeneratorAction);
 // runManager->SetUserAction(new RE04EventAction);
 // runManager->SetUserAction(new RE04TrackingAction);
 // runManager->SetUserAction(new RE04SteppingAction);
  
#ifdef G4VIS_USE
 // Visualization manager
 //
 G4VisManager* visManager = new G4VisExecutive;
 visManager->Initialize();
#endif
    
 // Initialize G4 kernel
 //
 runManager->Initialize();
  
 // Get the pointer to the User Interface manager
 //
 G4UImanager* pUImanager = G4UImanager::GetUIpointer();  

 if (argc==1)   // Define UI session for interactive mode
 {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      pUImanager->ApplyCommand("/control/execute vis.mac");
#endif
      ui->SessionStart();
      delete ui;
#endif
 }
 else           // Batch mode
 {
   G4String command = "/control/execute ";
   G4String fileName = argv[1];
   pUImanager->ApplyCommand(command+fileName);
 }

#ifdef G4VIS_USE
 delete visManager;
#endif
 delete runManager;

 return 0;
}
