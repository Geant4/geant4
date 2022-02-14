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
<<<<<<< HEAD
// $Id: examplePar01.cc 90093 2015-05-13 11:59:54Z gcosmo $
=======
//
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
//
//
// --------------------------------------------------------------
//      GEANT 4 - examplePar01
// --------------------------------------------------------------
// Comments
//
// Example of a main program making use of parameterisation
// i.e. "Fast Simulation"
//-------------------------------------------------------------------

//-------------------
// Generator action:
//-------------------
#include "Par01PrimaryGeneratorAction.hh"

//---------------------------
// Parameterisation manager:
//---------------------------
#include "G4GlobalFastSimulationManager.hh"

//------------
// Geometries:
//------------
#include "Par01DetectorConstruction.hh"
#include "Par01ParallelWorldForPion.hh"

//-----------------------------------
// Par01PhysicsList (makes use of the
// G4ParameterisationManagerProcess):
//-----------------------------------
#include "Par01PhysicsList.hh"

#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

// ----------------------------------------------------------------
// -- Action initialization (includes the primary generator action:
// ----------------------------------------------------------------
#include "Par01ActionInitialization.hh"

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
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  runManager->SetNumberOfThreads(4);

  // Detector/mass geometry and parallel geometry(ies):
  G4VUserDetectorConstruction* detector = new Par01DetectorConstruction();
  // -- Parallel geometry for pion parameterisation
  detector->RegisterParallelWorld(new Par01ParallelWorldForPion("pionGhostWorld"));
  // --  The name passed must be the same passed to the
  // -- G4FastSimulationManagerProcess attached to the pions
  runManager->SetUserInitialization(detector);
  
  // PhysicsList (including G4FastSimulationManagerProcess)
  G4VUserPhysicsList* physicsList = new Par01PhysicsList;
  runManager->SetUserInitialization(physicsList);

  //-------------------------------
  // UserAction classes
  //-------------------------------
  runManager->SetUserInitialization( new Par01ActionInitialization );

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

  if(argc==1)
  {
    //--------------------------
    // Define (G)UI
    //--------------------------
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete ui;
#endif
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
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
