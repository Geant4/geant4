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
// --------------------------------------------------------------
//      GEANT4 - OpticalEscape
//
// --------------------------------------------------------------
// Author: Tatiana Nikitina, CERN
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4ios.hh"

#include "AXPETDetectorConstruction.hh"
#include "AXPETPhysicsList.hh"
#include "AXPETPrimaryGeneratorAction.hh"
#include "AXPETRunAction.hh"
#include "AXPETEventAction.hh"
#include "AXPETSteppingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "Randomize.hh"


int main(int argc,char** argv)
{
  // Seed the random number generator manually
  //
  //G4long myseed = 345354;
  //CLHEP::HepRandom::setTheSeed(myseed);
 
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
   
  // Run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // UserInitialization classes - mandatory
  //
  AXPETDetectorConstruction* detector = new AXPETDetectorConstruction;
  runManager-> SetUserInitialization(detector);
  //
  AXPETPhysicsList* physics = new AXPETPhysicsList;
  runManager-> SetUserInitialization(physics);
  //
  AXPETPrimaryGeneratorAction* gen_action = new AXPETPrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);

#ifdef G4VIS_USE
  // visualization manager
  //
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  
  // UserAction classes
  //
  AXPETRunAction* run_action = new AXPETRunAction;
  runManager->SetUserAction(run_action);
  //
  AXPETEventAction* event_action = new AXPETEventAction();
  runManager->SetUserAction(event_action);
  //
  AXPETSteppingAction* step_action = new AXPETSteppingAction(detector);
  runManager->SetUserAction(step_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
    
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
   
  if (argc==1)   // Define UI session for interactive mode
    {
      G4UIsession* session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif    
      UI->ApplyCommand("/control/execute vis.mac"); 
      session->SessionStart();
      delete session;
   }
   
  else         // Batch mode
   {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
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

