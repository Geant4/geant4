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
//
//
// --------------------------------------------------------------
//      GEANT 4 - examplePar03
// --------------------------------------------------------------
// Comments
//
// Example of a main program demonstrating multiple energy deposition
// from fast simulation model.
//
//-------------------------------------------------------------------

#include "Par03DetectorConstruction.hh"
#include "Par03ActionInitialisation.hh"

#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4HadronicProcessStore.hh"
#include "G4EmParameters.hh"
#include "G4FastSimulationPhysics.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include <sstream>

int main(int argc, char** argv)
{
  // Macro name from arguments
  G4String batchMacroName;
  G4bool useInteractiveMode = true;
  G4String helpMsg(
    "Usage: " + G4String(argv[0]) +
    " [option(s)] \n No additional arguments triggers an interactive mode "
    "executing vis.mac macro. \n Options:\n\t-h\t\tdisplay this help "
    "message\n\t-m MACRO\ttriggers a batch mode executing MACRO\n");
  for(G4int i = 1; i < argc; ++i)
  {
    G4String argument(argv[i]);
    if(argument == "-h" || argument == "--help")
    {
      G4cout << helpMsg << G4endl;
      return 0;
    }
    else if(argument == "-m")
    {
      batchMacroName     = G4String(argv[i + 1]);
      useInteractiveMode = false;
      ++i;
    }
    else
    {
      G4Exception("main", "Unknown argument", FatalErrorInArgument,
                  ("Unknown argument passed to " + G4String(argv[0]) + " : " +
                   argument + "\n" + helpMsg)
                    .c_str());
    }
  }

  // Initialization of default Run manager
  auto* runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Detector geometry:
  auto detector = new Par03DetectorConstruction();
  runManager->SetUserInitialization(detector);

  // Physics list
  auto physicsList = new FTFP_BERT();
  // Add fast simulation physics
  auto fastSimulationPhysics = new G4FastSimulationPhysics();
  fastSimulationPhysics->BeVerbose();
  fastSimulationPhysics->ActivateFastSimulation("e-");
  fastSimulationPhysics->ActivateFastSimulation("e+");
  fastSimulationPhysics->ActivateFastSimulation("gamma");
  physicsList->RegisterPhysics(fastSimulationPhysics);
  // reduce verbosity of physics lists
  G4EmParameters::Instance()->SetVerbose(0);
  runManager->SetUserInitialization(physicsList);
  G4HadronicProcessStore::Instance()->SetVerbose(0);

  //-------------------------------
  // UserAction classes
  //-------------------------------
  runManager->SetUserInitialization(new Par03ActionInitialisation(detector));

  //----------------
  // Visualization:
  //----------------
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(useInteractiveMode)
  {
    auto  ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + batchMacroName);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}
