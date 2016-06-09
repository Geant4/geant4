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
// $Id: exampleN07.cc,v 1.12 2009-10-30 15:29:55 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN07 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "ExN07DetectorConstruction.hh"
#include "ExN07ParallelWorld.hh"
#include "ExN07PhysicsList.hh"
#include "ExN07PrimaryGeneratorAction.hh"
#include "ExN07RunAction.hh"
#include "ExN07SteppingVerbose.hh"

int main(int argc,char** argv)
{
 // Construct the stepping verbose class
 ExN07SteppingVerbose* verbosity = new ExN07SteppingVerbose;
 G4VSteppingVerbose::SetInstance(verbosity);

 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

 // Set mandatory initialization classes
 //
 G4VUserDetectorConstruction* detector = new ExN07DetectorConstruction;
 detector->RegisterParallelWorld(new ExN07ParallelWorld("ParallelScoringWorld"));
 runManager->SetUserInitialization(detector);
 //
 G4VUserPhysicsList* physics = new ExN07PhysicsList;
 runManager->SetUserInitialization(physics);
    
 // Set user action classes
 //
 G4VUserPrimaryGeneratorAction* gen_action = new ExN07PrimaryGeneratorAction;
 runManager->SetUserAction(gen_action);
 //
 G4UserRunAction* run_action = new ExN07RunAction;
 runManager->SetUserAction(run_action);
  
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
 G4UImanager* UImanager = G4UImanager::GetUIpointer();  

 if (argc==1)   // Define UI session for interactive mode
 {
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
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

#ifdef G4VIS_USE
 delete visManager;
#endif
 delete runManager;

 return 0;
}
