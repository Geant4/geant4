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
/// \file Par01/examplePar01.cc
/// \brief Main program of the Par01 example
//
//
//
// 
// --------------------------------------------------------------
//                     Geant4 - examplePar01
// --------------------------------------------------------------
// Comments
//
// Example of a main program making use of parameterisation
// i.e. "Fast Simulation"
//-------------------------------------------------------------------

#include "G4Types.hh"

//---------------
// -- Geometries:
//---------------
#include "Par01DetectorConstruction.hh"
#include "Par01ParallelWorldForPion.hh"

//---------------------------------------------------------------------
// -- Physics list, and tool to modify it and activate fast simulation:
//---------------------------------------------------------------------
#include "FTFP_BERT.hh"
#include "G4FastSimulationPhysics.hh"

#include "G4UImanager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

// ----------------------------------------------------------------
// -- Action initialization (includes the primary generator action:
// ----------------------------------------------------------------
#include "Par01ActionInitialization.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  //-------------------------------
  // Initialization of Run manager
  //-------------------------------
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4);
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|              Constructing MT run manager              |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
#else
  G4RunManager * runManager = new G4RunManager;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
  G4cout<<"|        Constructing sequential run manager            |"<<G4endl;
  G4cout<<"+-------------------------------------------------------+"<<G4endl;
#endif

  // -----------------------------------------------------
  // -- Detector/mass geometry and parallel geometry(ies):
  // -----------------------------------------------------
  // -- Mass geometry (ie : standard geometry):
  G4VUserDetectorConstruction* detector = new Par01DetectorConstruction();
  // -- Parallel geometry for pion parameterisation
  detector->RegisterParallelWorld(new Par01ParallelWorldForPion("pionGhostWorld"));
  // -- Passed to the run manager:
  runManager->SetUserInitialization(detector);

  // ----------------------------------------------
  // -- PhysicsList and fast simulation activation:
  // ----------------------------------------------
  // -- Create a physics list (note : FTFP_BERT is a G4VModularPhysicsList
  // -- which allows to use the subsequent G4FastSimulationPhysics tool to
  // -- activate the fast simulation):
  FTFP_BERT* physicsList = new FTFP_BERT;
  // -- Create helper tool, used to activate the fast simulation:
  G4FastSimulationPhysics* fastSimulationPhysics = new G4FastSimulationPhysics();
  fastSimulationPhysics->BeVerbose();
  // -- activation of fast simulation for particles having fast simulation models
  // -- attached in the mass geometry:
  fastSimulationPhysics->ActivateFastSimulation("e-");
  fastSimulationPhysics->ActivateFastSimulation("e+");
  fastSimulationPhysics->ActivateFastSimulation("gamma");
  // -- activation of fast simulation for particles having fast simulation models
  // -- attached in the parallel geometry:
  fastSimulationPhysics->ActivateFastSimulation("pi+","pionGhostWorld");
  fastSimulationPhysics->ActivateFastSimulation("pi-","pionGhostWorld");
  // -- Attach the fast simulation physics constructor to the physics list:
  physicsList->RegisterPhysics( fastSimulationPhysics );
  // -- Finally passes the physics list to the run manager:
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
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize ();

  if(ui)
  {
    //--------------------------
    // Define (G)UI
    //--------------------------
    ui->SessionStart();
    delete ui;
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

  delete visManager;
  delete runManager;

  return 0;
}
