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
// $Id: exampleN05.cc,v 1.17 2009-10-30 15:10:56 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN05
// --------------------------------------------------------------
// Comments
//
// Example of a main program making use of parameterisation
// i.e. "Fast Simulation"
//-------------------------------------------------------------------

//------------------------------
// Run, Generator etc.. actions:
//------------------------------
#include "ExN05RunAction.hh"
#include "ExN05PrimaryGeneratorAction.hh"
#include "ExN05EventAction.hh"
#include "ExN05SteppingAction.hh"

//---------------------------
// Parameterisation manager:
//---------------------------
#include "G4GlobalFastSimulationManager.hh"

//------------
// Geometries:
//------------
#include "ExN05DetectorConstruction.hh"
#include "ExN05ParallelWorldForPion.hh"

//-----------------------------------
// ExN05PhysicsList (makes use of the
// G4ParameterisationManagerProcess):
//-----------------------------------
#include "ExN05PhysicsList.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerKernel.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv)
{
  //-------------------------------
  // Initialization of Run manager
  //-------------------------------
  G4cout << "RunManager construction starting...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // Detector/mass geometry and parallel geometry(ies):
  G4VUserDetectorConstruction* detector = new ExN05DetectorConstruction();
  // -- Parallel geometry for pion parameterisation
  detector->RegisterParallelWorld(new ExN05ParallelWorldForPion("pionGhostWorld"));
  // --  The name passed must be the same passed to the
  // -- G4FastSimulationManagerProcess attached to the pions
  runManager->SetUserInitialization(detector);
  
  // PhysicsList (including G4FastSimulationManagerProcess)
  G4VUserPhysicsList* physicsList = new ExN05PhysicsList;
  runManager->SetUserInitialization(physicsList);

  //-------------------------------
  // UserAction classes
  //-------------------------------
  G4UserRunAction* run_action = new ExN05RunAction;
  runManager->SetUserAction(run_action);
  //
  G4VUserPrimaryGeneratorAction* gen_action = new ExN05PrimaryGeneratorAction;
  runManager->SetUserAction(gen_action);
  //
  G4UserEventAction* event_action = new ExN05EventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new ExN05SteppingAction;
  runManager->SetUserAction(stepping_action);

  // Initialize Run manager
  runManager->Initialize();

  //----------------
  // Visualization:
  //----------------
#ifdef G4VIS_USE
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize ();
#endif

  // Setup commands
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/Step/Verbose 0");
  UImanager->ApplyCommand("/tracking/Verbose 1");
  UImanager->ApplyCommand("/gun/particle e-");
  UImanager->ApplyCommand("/gun/energy 100 MeV");
  UImanager->ApplyCommand("/gun/direction 0 0 1");
  UImanager->ApplyCommand("/gun/position 0 0 0");
  UImanager->ApplyCommand("/gun/direction 0 .3 1.");

  if(argc==1)
  {
    //--------------------------
    // Define (G)UI
    //--------------------------
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
