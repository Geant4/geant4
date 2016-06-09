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
// $Id: pythia6_decayer.cc,v 1.1 2008/11/03 14:47:54 ivana Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
// ----------------------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

// ----------------------------------------------------------------------------

int main(int argc,char** argv)
{
  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
     
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new PhysicsList;
  runManager->SetUserInitialization(physics);
    
  // Set user action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = 
                          new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  RunAction* run_action = new RunAction;  
  runManager->SetUserAction(run_action);
  //
  EventAction* event_action = new EventAction(run_action);
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action =
                    new SteppingAction(detector, event_action);
  runManager->SetUserAction(stepping_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
  // Get the pointer to the User Interface manager
  //
  G4UImanager* UI = G4UImanager::GetUIpointer();      
  
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);    
    }
  else           // interactive mode : define visualization UI terminal
    {
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif

      G4UIsession* session = 0;
#if defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif

      UI->ApplyCommand("/control/execute vis.mac");
      session->SessionStart();
      delete session;
      
#ifdef G4VIS_USE
      delete visManager;
#endif                
    }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  delete runManager;

  return 0;
}

// ----------------------------------------------------------------------------
