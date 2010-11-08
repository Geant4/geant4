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
// $Id: RE03.cc,v 1.5 2010-11-08 20:14:41 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4ScoringManager.hh"

#include "RE03DetectorConstruction.hh"
#include "RE03PhysicsList.hh"
#include "RE03PrimaryGeneratorAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//====================================================================
// Un-comment this line for user defined score writer
//    #include "RE03UserScoreWriter.hh"
//====================================================================

int main(int argc,char** argv)
{
 // Construct the run manager
 //
 G4RunManager * runManager = new G4RunManager;

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
 G4VUserPhysicsList* physics = new RE03PhysicsList;
 runManager->SetUserInitialization(physics);
    
 // Set user action classes
 //
 G4VUserPrimaryGeneratorAction* gen_action = new RE03PrimaryGeneratorAction;
 runManager->SetUserAction(gen_action);
  
#ifdef G4VIS_USE
 // Visualization manager
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
   G4UIExecutive* ui = new G4UIExecutive(argc, argv);
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
