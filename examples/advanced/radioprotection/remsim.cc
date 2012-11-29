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

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "RemSimDetectorConstruction.hh"
#include "RemSimPhysicsList.hh"
#include "RemSimPrimaryGeneratorAction.hh"
#include "RemSimEventAction.hh"
#include "RemSimRunAction.hh"
#include "RemSimSteppingAction.hh"
#include "G4ScoringManager.hh"
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif
#include "QGSP_BIC_HP.hh"

int main(int argc,char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Access to the Scoring Manager pointer
  G4ScoringManager::GetScoringManager();

  //  HepRandom :: setTheSeed(0);
  // Set mandatory initialization classes
  
  // Geometry
  RemSimDetectorConstruction* detector = new RemSimDetectorConstruction();
  runManager -> SetUserInitialization(detector);
  
  // Physics
   runManager->SetUserInitialization(new RemSimPhysicsList);
  //runManager -> SetUserInitialization(new QGSP_BIC_HP());
  
  // Set mandatory user action class

  // Primary particles
  RemSimPrimaryGeneratorAction* primary = new RemSimPrimaryGeneratorAction();
  runManager -> SetUserAction(primary);
 
  // Set optional user action class
  RemSimEventAction* event = new RemSimEventAction();
  runManager -> SetUserAction(event);

  RemSimRunAction* run = new RemSimRunAction();
  runManager -> SetUserAction(run);

  runManager -> SetUserAction(new RemSimSteppingAction(primary));

#ifdef G4VIS_USE
   // Visualisation
   G4VisManager* visManager = new G4VisExecutive;
   visManager -> Initialize();
#endif
 
  // get the pointer to the UI manager and set verbosities
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
   
  if(argc == 1)
    // Define (G)UI terminal for interactive mode  
    { 
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute interactive.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  else
    // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager -> ApplyCommand(command+fileName);
    }


#ifdef G4VIS_USE
 delete visManager;
#endif

  // job termination
  delete runManager;
 
  return 0;
}

